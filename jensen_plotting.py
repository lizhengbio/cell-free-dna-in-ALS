import pickle
import matplotlib.pyplot as plt
import numpy as np

percent_meth = pickle.load(open("meth_dict.pkl", "rb"))

for k, v in percent_meth.items():
    percent_meth[k] = [round(x, 4) for x in v]


width = 1
plt.(list(range(101)), percent_meth["control"], width, color="blue")
plt.bar(list(range(101)), percent_meth["placenta"], width, color="red")
plt.bar(list(range(101)), percent_meth["buffy"], width, color="green")

plt.show()