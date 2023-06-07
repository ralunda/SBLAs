"""This piece of code is used to read the snapshots from HaloFinder outputs
to compile the snapshots file"""
import pandas as pd
import glob
import tqdm.notebook as tqdm

# directory of the HaloFinder results
path="/global/cfs/cdirs/agora/paper_CGM/Analysis/yt-4/ThinTimesteps/HaloFinder/rockstar_halos_all/ENZO/rockstar_halos/"

keep_columns=[
    "name", "galaxy_pos_x", "galaxy_pos_y", "galaxy_pos_z",
    "rho_max", "z_max", "M200[Msun]"]

files=glob.glob(f"{path}out*.list")
dfs = []
for file in tqdm.tqdm(files):
    f = open(file, "r")
    names = f.readline()[1:].split()
    aexpn = float(f.readline().split("=")[1])
    f.close()
    df=pd.read_csv(file, delim_whitespace=True, comment="#", names=names)
    num = int(file.replace(f"{path}out_", "").replace(".list", ""))
    df["name"] = f"RD{num:04d}"
    df["galaxy_pos_x"] = df["X"] / 0.7 * aexpn * 1000
    df["galaxy_pos_y"] = df["Y"] / 0.7 * aexpn * 1000
    df["galaxy_pos_z"] = df["Z"] / 0.7 * aexpn * 1000
    df["rho_max"] = df["Rvir"] / 0.7 * aexpn
    df["z_max"] = 1. / aexpn - 1.
    df["M200[Msun]"] = df["Mvir"] / 0.7
    dfs.append(df[keep_columns].copy())
final_df = pd.concat(dfs, ignore_index=True).sort_values("name")
final_df.to_csv(
    "/global/cfs/projectdirs/agora/ARRAKIHS/bin/SBLAs/Data/Rvir_cent_ENZO.txt.zip",
    sep=" ", index=False)
