from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

rows = []

with fits.open("J0559-1404_sequence.fits") as hdul:
    for hdu in hdul[1:]:
        if hdu.data is None:
            continue

        epoch = hdu.header.get("EPOCH")
        data = hdu.data

        for star_id, X, Y, ra, dec, mag, errmag in zip(
            data["ID_star"],
            data["X"],
            data["Y"],
            data["RA"],
            data["DEC"],
            data["mag"],
            data["errmag"]
        ):
            rows.append({
                "epoch": epoch,
                "star_id": int(star_id),
                "X": X,
                "Y": Y,
                "ra": (360/(2*np.pi)) * float(ra),
                "dec": (360/(2*np.pi)) * float(dec),
                "mag": mag,
                "errmag": errmag
            })

df_long = pd.DataFrame(rows)

# juntar ra y dec
#df_long[["ra","dec"]] = (
#    df_long.groupby("star_id")[["ra","dec"]]
#    .transform("mean")
#)

# calcular estadísticas por epoch
stats_epoch = df_long.groupby("epoch")["mag"].agg(
    median="median",
    std="std"
)

df_sigma = df_long.join(stats_epoch, on="epoch")


# sigma clip iterativo POSIBLE MEJORAR HACIENDOLO POR ESTRELLA EN VEZ DE POR EPOCH
#df_clip = df_long.copy()

#for _ in range(10):
#    stats = df_clip.groupby("epoch")["mag"].agg(median="median", std="std")
#    df_clip = df_clip.join(stats, on="epoch", rsuffix="_tmp")
#    df_clip = df_clip[
#        (df_clip["mag"] >= df_clip["median"] - 3 * df_clip["std"]) &
#        (df_clip["mag"] <= df_clip["median"] + 3 * df_clip["std"])
#    ]
#    df_clip = df_clip.drop(columns=["median", "std"])
#    print(len(df_clip))

# stats after clipping
#stats_epoch = df_clip.groupby("epoch")["mag"].agg(
#    median="median",
#    std="std"
#) 

#df_sigma = df_long.join(stats_epoch, on="epoch")

# extraer estrella
for star_id in df_sigma['star_id'].unique():
    #star_id = 2
    df_star1 = df_sigma[df_sigma["star_id"] == star_id].set_index("epoch")
    plt.errorbar(df_star1.index, df_star1['mag'] - df_star1['median'], yerr=df_star1['errmag'], fmt='--o')
    plt.xlabel("Epoch")
    plt.ylabel("Magnitude star - median")
    plt.title(f"Light curves")
    plt.gca().invert_yaxis()
plt.show()

star_id = 99
df_star1 = df_sigma[df_sigma["star_id"] == star_id].set_index("epoch")
plt.errorbar(df_star1.index, df_star1['mag'] - df_star1['median'], yerr=df_star1['errmag'], fmt='o')
plt.xlabel("Epoch")
plt.ylabel("Magnitude star - median")
plt.title(f"Light curve for star {star_id}")
plt.gca().invert_yaxis()
plt.show()

star_id = 102
df_star1 = df_sigma[df_sigma["star_id"] == star_id].set_index("epoch")
plt.errorbar(df_star1.index, df_star1['mag'] - df_star1['median'], yerr=df_star1['errmag'], fmt='o')
plt.xlabel("Epoch")
plt.ylabel("Magnitude star - median")
plt.title(f"Light curve for star {star_id}")
plt.gca().invert_yaxis()
plt.show()