# ------------------------------------------------------------
# Geometry-related utilitarians
# ------------------------------------------------------------

"""This is a modified version of:
   Protoro (Version 1.3)                                                     
   Generate 3D protonation states using Dimorphite_DL from a 3D sdf-structure provided as input_sdf                                                                 
   Author: Dr. Tony Froehlich 
   Protoro is using Dimorphite-DL (ver. 1.2.4), so please use the following citation:
   Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An open-source program 
   for enumerating the ionization states of drug-like small molecules. 
   J Cheminform 11:14. doi:10.1186/s13321-019-0336-9."""

import os
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, QED, Crippen, rdmolops, Descriptors, rdMolDescriptors
from rdkit.Geometry import Point3D
from rdkit.Chem.MolStandardize import rdMolStandardize
from util.dimorphite_dl_v124_TFmod import *
import numpy as np
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import io
import math

def print_header_protoro():
    # Citation information
    """Prints out header information."""
    print("\nProtoro is an extension of Dimorphite-DL (ver. 1.2.4), so please use the following citation:")
    print("Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An")
    print("open-source program for enumerating the ionization states of drug-like small")
    print("molecules. J Cheminform 11:14. doi:10.1186/s13321-019-0336-9.\n")

def load_mol(input_file):
    ext = os.path.splitext(input_file)[1]
    if ext == ".sdf":
        suppl = Chem.SDMolSupplier(input_file)
        return [mol for mol in suppl if mol is not None]
    elif ext == ".pdb":
        mol = Chem.MolFromPDBFile(input_file)
        return [mol] if mol else []
    else:
        raise ValueError("Unsupported file format. Only .sdf and .pdb files are supported.")

#write molecules to an SDF file, saving SMILES and name as properties
def write_sdf(mols, output_sdf, filename=None):
    with Chem.SDWriter(output_sdf) as writer:
        for mol in mols:
            if filename:
                mol.SetProp("_Name", filename)
            mol.SetProp("SMILES", Chem.MolToSmiles(mol))
            writer.SetForceV3000(False)                       #set to True to write SDF files in V3000 format instead of V2000
            writer.write(mol)

#add molecular properties such as MW, logP, TPSA, QED, etc.
def add_properties(mol, output_sdf):
    props = {'Formal_Charge': rdmolops.GetFormalCharge(mol),
             'Nr_Heavy_Atoms': mol.GetNumHeavyAtoms(),
             'MW': Descriptors.ExactMolWt(mol),
             'logP': Crippen.MolLogP(mol),
             'TPSA': rdMolDescriptors.CalcTPSA(mol),
             'QED': QED.qed(mol)}
    for prop_name, prop_value in props.items():
        mol.SetProp(prop_name,str(prop_value))
    writer = Chem.SDWriter(output_sdf)
    writer.SetForceV3000(False)                        #set to True to write SDF files in V3000 format instead of V2000
    writer.write(mol)
    writer.close()

def reorder_atoms(template_mol, target_mol, debug_log_path=None, tag=""):
    """Reorder atoms in target_mol to match template_mol.
       RDKit substructure matching can fail if charges/aromaticity/stereo flags differ
       after writing/reading SDF. We therefore try progressively normalized matches."""
    def _mol_summary(m):
        try:
            atoms = [a.GetSymbol() for a in m.GetAtoms()]
            arom = sum(1 for a in m.GetAtoms() if a.GetIsAromatic())
            charge = rdmolops.GetFormalCharge(m)
            bonds = []
            for b in m.GetBonds():
                bt = str(b.GetBondType())
                if b.GetIsAromatic():
                    bt = bt + "(arom)"
                bonds.append(bt)
            #bond histogram
            bh = {}
            for bt in bonds:
                bh[bt] = bh.get(bt, 0) + 1
            return {"n_atoms": m.GetNumAtoms(),
                    "n_bonds": m.GetNumBonds(),
                    "formal_charge": charge,
                    "n_arom_atoms": arom,
                    "atoms": atoms,
                    "bond_hist": bh,
                    "smiles": Chem.MolToSmiles(Chem.RemoveHs(m),isomericSmiles = True)}
        except Exception as e:
            return {"error": repr(e)}

    def _dbg(msg):
        if not debug_log_path:
            return
        try:
            with open(debug_log_path,"a",encoding = "utf-8") as fh:
                fh.write(msg.rstrip() + "\n")
        except Exception:
            pass

    #logging
    if debug_log_path:
        _dbg(f"\n=== reorder_atoms | {tag} ===")
        _dbg("TEMPLATE summary: " + str(_mol_summary(template_mol)))
        _dbg("TARGET   summary: " + str(_mol_summary(target_mol)))

    def _try_match(tmpl, targ, label):
        m = targ.GetSubstructMatch(tmpl)
        _dbg(f"Attempt {label}: match_len={len(m) if m else 0} "
             f"tmpl_atoms={tmpl.GetNumAtoms()} targ_atoms={targ.GetNumAtoms()}")
        return m if m and len(m) == tmpl.GetNumAtoms() else None

    #1) direct
    match = _try_match(template_mol,target_mol,"1_direct")
    if match:
        return Chem.RenumberAtoms(target_mol,list(match))
    #2) uncharged copies
    try:
        uncharger = rdMolStandardize.Uncharger()
        tmpl_u = uncharger.uncharge(Chem.Mol(template_mol))
        targ_u = uncharger.uncharge(Chem.Mol(target_mol))
        match_u = _try_match(tmpl_u,targ_u,"2_uncharged")
        if match_u:
            return Chem.RenumberAtoms(target_mol,list(match_u))
        #3) normalized: uncharged + noH + kekulize + remove stereo
        tmpl_n = Chem.RemoveHs(Chem.Mol(tmpl_u),implicitOnly = False)
        targ_n = Chem.RemoveHs(Chem.Mol(targ_u),implicitOnly = False)
        #kekulize to normalize aromatic representation (can throw; ignore if so)
        try:
            Chem.Kekulize(tmpl_n,clearAromaticFlags = True)
        except Exception:
            pass
        try:
            Chem.Kekulize(targ_n,clearAromaticFlags = True)
        except Exception:
            pass
        #remove stereochemistry that may differ after SDF roundtrip
        Chem.RemoveStereochemistry(tmpl_n)
        Chem.RemoveStereochemistry(targ_n)
        match_n = _try_match(tmpl_n, targ_n,"3_norm_kek_stereo_noH")
        if match_n:
            # match_n indices are in targ_n ordering which is identical to target_mol ordering
            return Chem.RenumberAtoms(target_mol,list(match_n))
    except Exception:
        pass
    #if we reach here, all matching attempts failed -> dump additional info
    if debug_log_path:
        try:
            base = os.path.splitext(debug_log_path)[0]
            Chem.MolToMolFile(template_mol, base + "_template.mol")
            Chem.MolToMolFile(target_mol, base + "_target.mol")
            _dbg(f"Wrote MOL dumps: {base}_template.mol and {base}_target.mol")
        except Exception as e:
            _dbg("Could not write MOL dumps: " + repr(e))
    raise ValueError("The target molecule does not match the template molecule (after normalization attempts).")

def embed_with_constraints(mol, input_file):
    template = load_mol(input_file)[0]
    template = force_remove_explicit_H_atoms(template) # remove leading explicit [H]    
    uncharger = rdMolStandardize.Uncharger()
    template_uncharged = uncharger.uncharge(template)
    template_noH = Chem.RemoveHs(template_uncharged,implicitOnly = False)
    template_noH_uncharged = uncharger.uncharge(template_noH)
    mol_H = Chem.AddHs(mol)
    GetFF = lambda m,confId = -1: AllChem.MMFFGetMoleculeForceField(m,AllChem.MMFFGetMoleculeProperties(m),confId = confId)
    mol_embed1 = AllChem.ConstrainedEmbed(mol_H,template_noH_uncharged,useTethers = True,randomseed = 0xf00d,getForceField = GetFF)
    return AllChem.ConstrainedEmbed(mol_embed1,template_uncharged,useTethers = True,randomseed = 0xf00d,getForceField = GetFF)

def calculate_mmff_energy(mol):
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        raise ValueError("MMFF parameters missing.")
    ff = AllChem.MMFFGetMoleculeForceField(mol,props)
    return ff.CalcEnergy()

def calculate_energy(mol):
    try:
        return calculate_mmff_energy(mol)
    except Exception as e:
        print(f"MMFF failed, falling back to UFF (no minimization): {e}")
        molH = Chem.AddHs(mol,addCoords = True)
        ff = AllChem.UFFGetMoleculeForceField(molH)
        return ff.CalcEnergy()

def protonate_mol2smiles(input_file,smi_file,pHmin,pHmax,pH,pHspread):
    input_mol = load_mol(input_file)
    protonated_mols = run_with_mol_list(input_mol,min_ph = pHmin,max_ph = pHmax)
    protonated_smiles = [Chem.MolToSmiles(mol) for mol in protonated_mols]
    print(f'Protonation States for pH {pH} ± {pHspread} as SMILES strings:\n{protonated_smiles}\n')
    with open(smi_file, 'w') as f:
        for smi in protonated_smiles:
            f.write(smi + '\n')
    return protonated_smiles

#draw a single molecule with a legend and optional property block below
def draw_molecule_with_props(mol, legend, properties = None, mol_size = (600, 600)):
    if mol is None:
        return Image.new("RGB",mol_size,"white")
    drawer = rdMolDraw2D.MolDraw2DCairo(mol_size[0],mol_size[1])
    opts = drawer.drawOptions()
    opts.bondLineWidth = 3
    opts.legendFontSize = 40
    rdMolDraw2D.PrepareAndDrawMolecule(drawer,mol,legend = legend)
    drawer.FinishDrawing()
    img_bytes = drawer.GetDrawingText()
    mol_img = Image.open(io.BytesIO(img_bytes)).convert("RGB")
    #draw property block below if any properties given
    if properties:
        line_height = 35                        #increase value if font size was increased
        text_block_height = len(properties)*line_height + 20
        total_height = mol_size[1] + text_block_height
        final_img = Image.new("RGB",(mol_size[0],total_height),"white")
        final_img.paste(mol_img,(0,0))
        draw = ImageDraw.Draw(final_img)
        try:
            font_regular = ImageFont.truetype("arial.ttf",30)              #change value to adjust font size
        except:
            font_regular = ImageFont.load_default()
        y_text = mol_size[1] + 10
        for key, val in properties:
            draw.text((10,y_text),f"{key}: {val}",fill = "black",font = font_regular)
            y_text += line_height
        return final_img
    else:
        return mol_img

#draws molecules in a grid with legends and (optional) properties below each molecule
def draw_molecule_grid(mols, legends, prop_lists = None, molsPerRow = 4, mol_size = (600,600), out_path = None):
    num_mols = len(mols)
    n_rows = math.ceil(num_mols/molsPerRow)
    #draw all molecule images
    mol_imgs = []
    for i in range(num_mols):
        props = prop_lists[i] if prop_lists else None
        img = draw_molecule_with_props(mols[i],legends[i],props,mol_size)
        mol_imgs.append(img)
    #determine max height per row
    row_heights = []
    for r in range(n_rows):
        imgs_in_row = mol_imgs[r*molsPerRow:(r + 1)*molsPerRow]
        row_heights.append(max(img.height for img in imgs_in_row))
    #build grid
    grid_width = molsPerRow*mol_size[0]
    grid_height = sum(row_heights)
    grid_img = Image.new("RGB",(grid_width,grid_height),"white")
    y = 0
    for r in range(n_rows):
        imgs_in_row = mol_imgs[r*molsPerRow:(r + 1)*molsPerRow]
        x = 0
        row_height = row_heights[r]
        for img in imgs_in_row:
            grid_img.paste(img,(x,y))
            x += mol_size[0]
        y += row_height
    if out_path:
        grid_img.save(out_path)
        print(f"High-resolution 2D representation of protonation states saved to {out_path}")
    return grid_img

def force_remove_explicit_H_atoms(mol):
    if mol is None:
        return None
    rw = Chem.RWMol(mol)
    h_ids = [a.GetIdx() for a in rw.GetAtoms() if a.GetSymbol() == "H"]
    for idx in sorted(h_ids,reverse = True):
        rw.RemoveAtom(idx)
    m2 = rw.GetMol()
    m2.UpdatePropertyCache(strict = False)
    return m2

def RDKitProtonation(input_file, pH, image_mode, pH_spread = 0.5, logger = None, warning_log_path = None):
    """ORIGINAL FEATURES PRESERVED:
         - generates *_prot_*_unranked.sdf intermediates
         - ranks by MMFF energy, keeps template state as rank 1 when present
         - renames ranked outputs to:
               rank 1:  <base_name>.sdf
               rank k:  <base_name>_prot_<k>.sdf
         - generates images exactly like original (highres/lowres/none + filenames/sizes)
         - removes intermediate unranked SDF files (raw_files cleanup)
       NEW FEATURES ADDED:
         - Attempt 1: standard ConstrainedEmbed (original)
         - If that fails for a state: ETKDG fallback (per-state)
         - logger.warning + RMSD + reliable/moderate/less reliable label
         - append same message to <OUTPUT>/WARNING_ligand_input.log (only when fallback used)"""
    import datetime
    from rdkit.Chem import rdMolAlign

    def _warn(msg: str) -> None:
        #logger.warning
        try:
            if logger is not None:
                logger.warning(msg)
        except Exception:
            pass
        #append to WARNING file only if caller provided a path
        if warning_log_path:
            try:
                with open(warning_log_path,"a",encoding = "utf-8") as fh:
                    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    fh.write(f"[{ts}] {msg}\n")
            except Exception:
                pass

    def _rmsd_label(rmsd: float) -> str:
        #tiers you requested
        if rmsd <= 0.50:
            return "reliable"
        if rmsd <= 1.00:
            return "moderate"
        return "less reliable"

    def _fallback_embed_etkdg(mol_noH, template_mol):
        """Fallback: AddHs -> ETKDG -> MMFF optimize -> return a NO-H molecule
           (to match the original Protoro pipeline expectations).
           Also compute RMSD on uncharged NO-H copies for a stable match."""
        uncharger = rdMolStandardize.Uncharger()
        #template reference for RMSD (uncharged, noH)
        template_unch = uncharger.uncharge(template_mol)
        template_noH = Chem.RemoveHs(template_unch,implicitOnly = False)
        #build H-containing mol for embedding/optimization
        mH = Chem.AddHs(mol_noH)
        params = AllChem.ETKDGv3()
        params.randomSeed = 1
        if hasattr(params,"maxAttempts"):
            try:
                params.maxAttempts = 1000
            except Exception:
                pass
        status = AllChem.EmbedMolecule(mH,params)
        if status != 0:
            raise ValueError("ETKDG embedding failed (status != 0).")
        AllChem.MMFFOptimizeMolecule(mH,mmffVariant = "MMFF94s",maxIters = 500)
        #return NO-H molecule, because the pipeline expects to reorder atoms against a NO-H template
        m_noH = Chem.RemoveHs(mH,implicitOnly = False)
        #RMSD: use uncharged NO-H copies to avoid match failures due to formal charge differences
        try:
            m_noH_unch = uncharger.uncharge(Chem.Mol(m_noH))
            rmsd = rdMolAlign.AlignMol(m_noH_unch,template_noH)
            label = _rmsd_label(rmsd)
        except Exception:
            rmsd = float("nan")
            label = "unknown"
        return m_noH,rmsd,label
    #original body
    max_prots = 15
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_dir = input_file.replace(".sdf","").replace(".pdb","").replace(base_name,"")
    pHmin = pH - pH_spread
    pHmax = pH + pH_spread
    smi_file = os.path.join(output_dir,f"{base_name}.smi")
    smiles_list = protonate_mol2smiles(input_file,smi_file,pHmin,pHmax,pH,pH_spread)
    raw_files = []
    ranked_data = []
    #load template once (used for fallback RMSD reference)
    template_input = load_mol(input_file)[0]
    for i,smiles in enumerate(smiles_list):
        state_idx = i + 1
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("MolFromSmiles returned None")
            #PATCH: remove standalone explicit [H] atoms that Dimorphite may put in SMILES
            #example problematic SMILES patterns: "[H]/[NH+]=..." -> explicit H atom + NH implies another H
            mol = force_remove_explicit_H_atoms(mol)
            mol.UpdatePropertyCache(strict = False)
            # Attempt 1: original constrained embed
            try:
                embedded = embed_with_constraints(mol,input_file)
            except Exception as e_ce:
                #attempt 2: fallback ETKDG
                embedded,rmsd,label = _fallback_embed_etkdg(mol,template_input)
                _warn(f"ConstrainedEmbed failed. Used ETKDG fallback for state {state_idx}. "
                      f"RMSD_to_template(heavy)={rmsd:.3f} Å -> {label}. "
                      f"SMILES = {smiles} | error = {repr(e_ce)}")
            #IMPORTANT: downstream reorder_atoms assumes NO-H
            embedded = Chem.RemoveHs(embedded,implicitOnly = False)
            output_file = os.path.join(output_dir,f"{base_name}_prot_{state_idx}_unranked.sdf")
            write_sdf([embedded],output_file)
            raw_files.append(output_file)
            template_raw = load_mol(input_file)[0]
            if template_raw is None:
                raise ValueError("Template ligand could not be loaded.")
            #FORCE remove explicit 'H' atoms from the template for matching
            #this fixes cases like '[H]/N=C...' where RDKit keeps an explicit H atom.
            template_match = force_remove_explicit_H_atoms(template_raw)
            if template_match is None:
                raise ValueError("force_remove_explicit_H_atoms(template_raw) returned None")
            if template_match.GetNumConformers() == 0:
                raise ValueError("template_match has no conformer; cannot extract 3D coords.")
            #for matching: uncharge but KEEP ATOM ORDER (no SMILES round-trip!)
            uncharger = rdMolStandardize.Uncharger()
            template_for_match = uncharger.uncharge(Chem.Mol(template_match))
            #coordinates: extract from template_match in the SAME atom order as template_for_match
            conf_raw = template_match.GetConformer()
            coords = [(conf_raw.GetAtomPosition(i).x,
                       conf_raw.GetAtomPosition(i).y,
                       conf_raw.GetAtomPosition(i).z)
                      for i in range(template_match.GetNumAtoms())]
            #target: ensure explicit 'H' removed too
            mol_conf = force_remove_explicit_H_atoms(Chem.Mol(embedded))
            #reorder so target atom order matches template_for_match atom order
            reordered = reorder_atoms(template_for_match,mol_conf,debug_log_path = None,tag = f"{base_name} state {state_idx}")
            #now coords and reordered indices are consistent
            conf = reordered.GetConformer()
            if reordered.GetNumAtoms() != len(coords):
                raise ValueError(f"Coord mismatch: reordered has {reordered.GetNumAtoms()} atoms, "
                                 f"coords has {len(coords)}.")
            for j,(x,y,z) in enumerate(coords):
                conf.SetAtomPosition(j,Point3D(x,y,z))
            with_H = Chem.AddHs(reordered,addCoords = True)
            write_sdf([with_H],output_file)
            energy = calculate_energy(with_H)
            charge = rdmolops.GetFormalCharge(with_H)
            ranked_data.append((with_H,output_file,energy,charge,smiles))
        except Exception as e:
            print(f"Skipping prot {state_idx} due to error: {e}")
            #if we already wrote output_file, keep a copy for debugging (only when something goes wrong)
            try:
                if 'output_file' in locals() and output_file and os.path.exists(output_file):
                    debug_copy = os.path.join(output_dir,f"FAILED_{base_name}_prot_{state_idx}_unranked.sdf")
                    try:
                        os.replace(output_file,debug_copy)
                    except Exception:
                        #if replace fails (cross-FS), try copy
                        import shutil
                        shutil.copyfile(output_file,debug_copy)
            except Exception:
                pass
            continue
    template_mol = load_mol(input_file)[0]
    template_smiles = Chem.MolToSmiles(Chem.RemoveHs(template_mol))
    original_idx = None
    for i,(mol,file,energy,charge,smiles) in enumerate(ranked_data):
        smiles2 = Chem.MolToSmiles(Chem.RemoveHs(mol))
        if smiles2 == template_smiles:
            original_idx = i
            break
    original_prot = None
    if original_idx is not None:
        original_prot = ranked_data.pop(original_idx)
    sorted_prots = sorted(ranked_data,key = lambda x: x[2])[:max_prots - 1]
    if original_prot is not None:
        sorted_prots = [original_prot] + sorted_prots
    molecules,legends = [],[]
    listSDFfiles = []
    listCharges = []
    listSMILES = []
    for rank, (mol,file,energy,charge,smiles) in enumerate(sorted_prots,start = 1):
        if rank == 1 and original_prot is not None:
            sorted_filename = os.path.join(output_dir,f"{base_name}.sdf")
            prot_name = f"{base_name}"
        else:
            sorted_filename = os.path.join(output_dir,f"{base_name}_prot_{rank}.sdf")
            prot_name = f"{base_name}_prot_{rank}"
        os.rename(file,sorted_filename)
        if os.path.exists(file):
            os.remove(file)
        mol.SetProp("_Name",prot_name)
        mol_with_h = Chem.AddHs(mol,addCoords = True)
        write_sdf([mol_with_h],sorted_filename)
        mol_2d = Chem.RemoveHs(mol)
        smiles2 = Chem.MolToSmiles(mol_2d)
        AllChem.Compute2DCoords(mol_2d)
        molecules.append(mol_2d)
        legends.append(f"{prot_name}\nEnergy: {energy:.2f}\nCharge: {charge}")
        listSDFfiles.append(sorted_filename)
        listCharges.append(charge)
        listSMILES.append(smiles2)
    #image generation options (EXACT original behaviour)
    if image_mode == "highres":
        draw_molecule_grid(molecules,legends,
                           prop_lists = None,
                           molsPerRow = 4,
                           mol_size = (800,800),
                           out_path = os.path.join(output_dir,f"{base_name}_prots_highres.png"))
    elif image_mode == "lowres":
        img = Draw.MolsToGridImage(molecules,legends = legends,molsPerRow = 4,subImgSize = (300,300))
        image_file = os.path.join(output_dir,f"{base_name}_prots.png")
        img.save(image_file)
    elif image_mode == "none":
        pass
    #remove intermediate unranked SDF files (EXACT original behaviour)
    for raw_file in raw_files:
        if os.path.exists(raw_file):
            os.remove(raw_file)
    return listSDFfiles,listCharges,listSMILES
