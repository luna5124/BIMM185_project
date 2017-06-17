import matplotlib.pyplot as plt

distances = []
rhos = []
unique_rhos = []


with open ('dis_rho_new.out','r') as file:
	for line in file:
		line = line.strip().split(' ')
		distances.append(line[0])
		rhos.append(line[1])
		if line[1] not in unique_rhos:
			unique_rhos.append(line[1])

unique_rhos = sorted(unique_rhos,reverse=True)

ranks = []
for i in range(len(distances)):
	rho_index = unique_rhos.index(rhos[i])+1
	ranks.append(rho_index)
	print(distances[i], rho_index)

plt.plot(distances, ranks,'bo')
plt.show()





