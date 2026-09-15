#start from Miniconda
FROM condaforge/miniforge3:latest

#set working directory
WORKDIR /app

#copy repo into the image
COPY . .

#create conda environment with needed packages
RUN conda create -n myenv -c conda-forge -y \
    python=3.10 \
    rdkit \
    gemmi \
    mrcfile \
 && conda clean -a -y

#activate environment automatically
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

#install C++ compiler
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ \
 && rm -rf /var/lib/apt/lists/*

#compile .cpp:
RUN echo "Compiling gfn2xtb_in_pocket_opt..."
RUN g++ -O3 gfn2xtb_in_pocket_opt.cpp -o gfn2xtb_in_pocket_opt.exe
RUN echo "Compiling align..."
RUN g++ -O3 align.cpp -o align.exe
RUN echo "Compiling CCD_Protonate..."
RUN g++ -O3 CCD_Protonate.cpp -o CCD_Protonate.exe
RUN echo "Compiling EDA_GFN2_3mol_solv..."
RUN g++ -O3 EDA_GFN2_3mol_solv.cpp -o EDA_GFN2_3mol_solv.exe
RUN echo "Compiling Protein_Protonate..."
RUN g++ -O3 Protein_Protonate.cpp -o Protein_Protonate.exe
RUN echo "Compiling rmsd..."
RUN g++ -O3 rmsd.cpp -o rmsd.exe
RUN echo "Compiling VSEPR_Protonate..."
RUN g++ -O3 VSEPR_Protonate.cpp -o VSEPR_Protonate.exe
RUN echo "Compiling gfn2xtb_density..."
RUN g++ -O3 gfn2xtb_density.cpp -o gfn2xtb_density.exe
RUN echo "Compiling SDFMOL2_Protonate..."
RUN g++ -O3 SDFMOL2_Protonate.cpp -o SDFMOL2_Protonate.exe
RUN echo "Compiling gfn2xtb_in_pocket_general_opt..."
RUN g++ -O3 gfn2xtb_in_pocket_general_opt.cpp -o gfn2xtb_in_pocket_general_opt.exe
RUN echo "All C++ files compiled successfully."
#default command
CMD ["conda", "run", "--no-capture-output", "-n", "myenv", "python", "in_pocket.py"]
