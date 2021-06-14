from nxcode import readx
import matplotlib.pyplot as plt
import networkx as nx

f = plt.figure(figsize=(35,25))

plt.subplot(2, 3, 1)
G = readx("../../magna_UPGMA/Data/a.txt")
d = dict(G.degree)
nx.draw_networkx(G, with_labels = True, node_color="#d9ffb3", edgecolors='black', node_size=[v * 500 for v in d.values()])
plt.title('Wuttagoonaspis (a) \n', fontsize=30)
plt.axis('off')

plt.subplot(2, 3, 2)
G = readx("../../magna_UPGMA/Data/b.txt")
d = dict(G.degree)
nx.draw_networkx(G, with_labels = True, node_color="#d9ffb3", edgecolors='black', node_size=[v * 500 for v in d.values()])
plt.title('Romundina (b) \n', fontsize=30)
plt.axis('off')

plt.subplot(2, 3, 3)
G = readx("../../magna_UPGMA/Data/c.txt")
d = dict(G.degree)
nx.draw_networkx(G, with_labels = True, node_color="#d9ffb3", edgecolors='black', node_size=[v * 500 for v in d.values()])
plt.title('Brindabellaspis (c) \n', fontsize=30)
plt.axis('off')

plt.subplot(2, 3, 4)
G = readx("../../magna_UPGMA/Data/d.txt")
d = dict(G.degree)
nx.draw_networkx(G, with_labels = True, node_color="#d9ffb3", edgecolors='black', node_size=[v * 500 for v in d.values()])
plt.title('Eurycaraspis (d) \n', fontsize=30)
plt.axis('off')

plt.subplot(2, 3, 5)
G = readx("../../magna_UPGMA/Data/e.txt")
d = dict(G.degree)
nx.draw_networkx(G, with_labels = True, node_color="#d9ffb3", edgecolors='black', node_size=[v * 500 for v in d.values()])
plt.title('Dicksonosteus (e) \n', fontsize=30)
plt.axis('off')

plt.subplot(2, 3, 6)
G = readx("../../magna_UPGMA/Data/f.txt")
d = dict(G.degree)
nx.draw_networkx(G, with_labels = True, node_color="#d9ffb3", edgecolors='black', node_size=[v * 500 for v in d.values()])
plt.title('Entelognathus (f) \n', fontsize=30)
plt.axis('off')




plt.savefig("Placoderms.png")
plt.close()