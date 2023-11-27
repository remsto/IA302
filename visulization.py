import matplotlib.pyplot as plt
import numpy as np

def create_CO2():
    atoms = np.array([[-1, 0, 0],
                      [0, 0, 0],
                      [1, 0, 0]])
    links = np.array([[0, 1, 0],
                      [1, 0, 1],
                      [0, 1, 0]])
    
    return atoms, links

def create_CO2_xyz():
    str = ""
    str += "3\n"
    str += "Carbon Dioxyde Molecule"
    str += "C  0.000  0.000  0.000"
    str += "O  -1.000  0.000  0.000"
    str += "O  1.000  0.000  0.000"

def parse_results(filename):
    """
    Open a .txt file of format
    ([-0.05000000000000001, 0.05000000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [0.95, 1.050000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [-0.03572656471497753, 0.0140459505491084] ; [-0.9208919808319979, -0.8629367518976043] ; [-0.05000000000000001, 0.05000000000000001] ; [0.000908437500000303, 0.09580202578125036] ; [0.4695187500000004, 0.5644123382812506] ; [0.7325848641796881, 0.7963744429687508])
    ...
    ([-0.05000000000000001, 0.05000000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [0.95, 1.050000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [0.01404595054910839, 0.07487902476076897] ; [-0.9208919808319979, -0.8629367518976043] ; [-0.05000000000000001, 0.05000000000000001] ; [0.000908437500000303, 0.09580202578125036] ; [0.4695187500000004, 0.5644123382812506] ; [0.7325848641796881, 0.7963744429687508])
    
    and return a list of strings representing the individual solutions to the problem
    """
    with open(filename, "r") as f:
        file = f.read()
    i = file.find("Uncertain intervals:\n")
    file = file[i + 21:].split("\n")
    return file

def get_barycentres(sol):
    """
    Open a solution string of format
    ([-0.05000000000000001, 0.05000000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [0.95, 1.050000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [-0.05000000000000001, 0.05000000000000001] ; [-0.03572656471497753, 0.0140459505491084] ; [-0.9208919808319979, -0.8629367518976043] ; [-0.05000000000000001, 0.05000000000000001] ; [0.000908437500000303, 0.09580202578125036] ; [0.4695187500000004, 0.5644123382812506] ; [0.7325848641796881, 0.7963744429687508])
    
    and return the list of the center of each interval
    """
    sol = sol[1:-1]
    interstrings = sol.split(" ; ")
    inters = [s[1:-1].split(", ") for s in interstrings]
    bars = [float(inter[0]) + (float(inter[1]) - float(inter[0]))/2 for inter in inters]
    return bars

def get_average_solution(sols):
    """
    Open the file of solutions
    and return ONE solution corresponding to the barycenter of all solutions
    """
    solutions = []
    for sol in sols:
        if sol:
            bars = np.array(get_barycentres(sol))
            solutions.append(bars)
    solutions = np.asarray(solutions)
    average = np.mean(solutions, axis=0)
    print(average)
    print(solutions[1000])
    return average.tolist()
        

def barycentres_to_xyz(bars):
    """
    Open a list of coordinates of format
    [x0, y0, z0, x1, y1, ..., xn, yn, zn]
    
    and return a string corresponding the .xyz file format
    """
    xyz = ""
    for k in range(0, len(bars), 3):
        xyz += ("C " + str(bars[k]) + " " + str(bars[k+1]) + " " + str(bars[k+2]) + "\n")
    return xyz



def pltvisualize(atoms, links):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Drawing atoms
    ax.scatter(atoms[:, 0], atoms[:, 1], atoms[:, 2])
    # Drawing links
    for i in range(len(atoms)):
        for j in range(i, len(atoms)):
            if links[i, j]:
                ax.plot(xs=[atoms[i][0], atoms[j][0]], ys=[atoms[i][1], atoms[j][1]], zs=[atoms[i][2], atoms[j][2]])

    plt.show()
    
if __name__=="__main__":
    if False:
        examples = parse_results("results.txt")
        bars = get_barycentres(examples[4])
        xyz = barycentres_to_xyz(bars)
        with open("pyramid.xyz", "w") as f:
            f.write(xyz)
    if False:
        examples = parse_results("results3.txt")
        for i, e in enumerate(examples):
            bars = get_barycentres(e)
            xyz = barycentres_to_xyz(bars)
            with open(f"pyramids3/pyramid{i}.xyz", "w") as f:
                f.write(xyz)
                
    if True:
        examples = parse_results("results3.txt")
        print(examples)
        average = get_average_solution(examples)
        xyz = barycentres_to_xyz(average)
        with open(f"pyramids3/pyramid_average.xyz", "w") as f:
            f.write(xyz)
        
