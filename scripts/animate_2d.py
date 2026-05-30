import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("output_0000.csv")

fig, ax = plt.subplots()

sc = ax.scatter(df["x"], df["y"], c=df["rho"], s=10)
fig.colorbar(sc, ax=ax, label="density")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")
plt.show()
