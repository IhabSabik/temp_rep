import sys

from skimage import io
from Functions import *

start_time = time()

voxels = io.imread('labeled_pores.tif')

dist_path = 'computed_distances.txt'

f = open(dist_path, 'a')
f.close()
log_file = open('log_file.txt', 'a+')
log_file_detailed = open('log_file_detailed.txt', 'a+')

with open('Centroids', 'r') as f:
    Centroids = [[float(num) for num in line.split(',')] for line in f]
num_pores = len(Centroids)

coordinates = [[] for item in range(num_pores)]

with open('list_boundary_alpha_shape.txt', 'r') as f:
    list_boundary = [int(elem) - 1 for elem in f.read().split(',') if elem]

with open('list_boundary_voxelized.txt', 'r') as f:
    list_boundary_vox = [int(elem) - 1 for elem in f.read().split(',') if elem]

with open('Diameters.txt', 'r') as f:
    Radii = [float(elem) / 2 for elem in f.read().split('\n') if elem]
ignore = [vertex for vertex in range(num_pores) if Radii[vertex] < (3 * 450 / (4 * math.pi)) ** (1 / 3)]

list_interior = set_diff(set_diff(range(num_pores), ignore), list_boundary)

center = sum(np.array(Centroids)) / len(Centroids)
dists_to_center = [np.linalg.norm(item - center) for item in Centroids]
indices_sorted = list(np.argsort(np.array(dists_to_center)))
lst_1 = [int(item) for item in indices_sorted if item in list_interior]
lst_2 = [int(item) for item in indices_sorted if item in list_boundary]
ordered_pores = lst_1 + lst_2

max_num_open_edges = 25
index_inf = math.inf
threshold = 1

Register_Added = []
Register_Removed = []

Matrix_out = [[] for item in Centroids]
avg_deg = 0

log_file_detailed.write("Compute initial adjacencies...")
log_file_detailed.flush()
log_file.write("Compute initial adjacencies...")
log_file.flush()
while avg_deg < 6:
    threshold += 1
    for vertex in set_diff(range(num_pores), ignore):
        for vertex_sub in set_diff(range(vertex + 1, num_pores), ignore):
            if vertex_sub in Matrix_out[vertex]:
                continue
            print("\rCompute initial adjacencies (tolerance = {}):".format(threshold), vertex, "out of", num_pores - 2,
                  " | ", vertex_sub, "out of", num_pores - 1, end="")
            if norm(np.array(Centroids[vertex]) - np.array((Centroids[vertex_sub]))) > Radii[vertex] + Radii[
                vertex_sub] + 10 * threshold:
                continue
            dist = distance_vox(vertex, vertex_sub, Centroids, coordinates, voxels, dist_path)
            if dist < threshold:
                Matrix_out = add_edge(vertex, vertex_sub, Matrix_out)
    avg_deg = 0
    for vertex in list_interior:
        avg_deg += len(Matrix_out[vertex])
    avg_deg /= len(list_interior)
    print("\nThreshold = {}  |  average face degree = {}".format(threshold, round(avg_deg, 3)))
    log_file.write("\nThreshold = {}  |  average face degree = {}".format(threshold, round(avg_deg, 3)))
    log_file.flush()
    log_file_detailed.write("\nThreshold = {}  |  average face degree = {}".format(threshold, round(avg_deg, 3)))
    log_file_detailed.flush()
print("")

log_file_detailed.write("\n\nCheck for isolated vertices...")
log_file_detailed.flush()
log_file.write("\n\nCheck for isolated vertices...")
log_file.flush()
isolated = [item for item in set_diff(range(num_pores), ignore) if not Matrix_out[item]]
print("Initially isolated:", len(isolated))
log_file_detailed.write("\nInitially isolated: {}".format(len(isolated)))
log_file_detailed.flush()
log_file.write("\nInitially isolated: {}".format(len(isolated)))
log_file.flush()
count = 0
for vertex in isolated:
    count += 1
    print("\rAdd connections to isolated vertices: processing {} out of {}".format(count, len(isolated)), end="")
    distances, vts = [], []
    for vertex_sub in set_diff(range(num_pores), ignore):
        if vertex_sub == vertex:
            continue
        if norm(np.array(Centroids[vertex]) - np.array((Centroids[vertex_sub]))) > Radii[vertex] + Radii[
            vertex_sub] + 15 * threshold:
            continue
        dist = distance_vox(vertex, vertex_sub, Centroids, coordinates, voxels, dist_path)
        distances.append(dist)
        vts.append(vertex_sub)
    m = min(distances)
    vertex_sub = vts[distances.index(m)]
    Matrix_out = add_edge(vertex, vertex_sub, Matrix_out)
print("")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
avg_deg /= len(list_interior)
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg, 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg, 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

# Add connections to endpoints of edges with low degree
start_time = time()
print("\nAdd connections to endpoints of edges with low degree...")
log_file_detailed.write("\n\nAdd connections to endpoints of edges with low degree...")
log_file_detailed.flush()
log_file.write("\n\nAdd connections to endpoints of edges with low degree...")
log_file.flush()
for u in range(num_pores):
    for v in [item for item in Matrix_out[u] if item > u]:
        if sublist([u, v], list_boundary):
            continue
        while len(mutual_neighbors([u, v], Matrix_out)) < 3:
            pairs = []
            for w in set_diff(Matrix_out[u], [v] + Matrix_out[v]):
                pairs.append([v, w])
            for w in set_diff(Matrix_out[v], [u] + Matrix_out[u]):
                pairs.append([u, w])
            if not pairs:
                break
            distances = [distance_vox(item[0], item[1], Centroids, coordinates, voxels, dist_path) for item in pairs]
            Ind = distances.index(min(distances))
            if distances[Ind] > 2 * threshold:
                break
            Register_Added.append(pairs[Ind])
            Matrix_out = add_edge(pairs[Ind][0], pairs[Ind][1], Matrix_out)
            print("Added:", pairs[Ind], distances[Ind])
            log_file_detailed.write("\nAdded: {} with distance {}".format(pairs[Ind], distances[Ind]))
            log_file_detailed.flush()

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

start_time = time()
print("\n\nCompute the square holes in vertex links...")
log_file_detailed.write("\n\nCompute the square holes in vertex links...")
log_file_detailed.flush()
log_file.write("\n\nCompute the square holes in vertex links...")
log_file.flush()
computed_holes = [[] for item in range(num_pores)]
deg_3_edges_in_cross = [[] for item in range(num_pores)]
computed_open_edges = [[] for item in range(num_pores)]
for vertex in range(num_pores):
    print("\rCompute the square holes in vertex links:", vertex, "out of", num_pores - 1, end="")
    Holes, edges_in_cross, open_edges = holes(vertex, Matrix_out, list_boundary, Centroids, True, holes_degree=4)
    computed_holes[vertex] = Holes
    deg_3_edges_in_cross[vertex] = edges_in_cross
    computed_open_edges[vertex] = open_edges
print("")

log_file_detailed.write("\n\nClose the square holes in vertex links...")
log_file_detailed.flush()
log_file.write("\n\nClose the square holes in vertex links...")
log_file.flush()
iteration = 0
while True:
    out = True
    iteration += 1
    index = -1
    while index < len(ordered_pores) - 1:
        index += 1
        print("\r", "Close the square holes:", index, "out of", num_pores - 1, "iteration", iteration, end="")
        vertex = ordered_pores[index]
        if len(computed_open_edges[vertex]) < max_num_open_edges:  # only close squared holes in complex links
            continue
        Holes = loads(dumps(computed_holes[vertex]))
        for i in range(len(Holes)):
            if len(Holes[i]) != 4:
                continue
            priority = []
            for j in range(len(Holes[i])):
                if is_line([Holes[i][item] for item in [j, Next(j, 1, len(Holes[i]))]],
                          deg_3_edges_in_cross[vertex]) and is_line(
                    [Holes[i][item] for item in [j, Previous(j, 1, len(Holes[i]))]], deg_3_edges_in_cross[vertex]):
                    priority.append(Holes[i][j])
            for j in range(len(Holes[i]) - 2):  # extract priority vertices
                if Holes[i][j] == Holes[i][j + 2]:
                    priority = union(priority, Holes[i][j + 1])
            indices = []
            for j in combinations(range(len(Holes[i])), 2):
                j = list(j)
                edge = [Holes[i][item] for item in j]
                if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                    continue
                if priority and not intersect(edge, priority) and len(Holes[i]) != 4:
                    continue
                indices.append(j)
            if not indices:
                continue

            distances, outer_vertices, other_edges = [], [], []
            for j in indices:
                edge = [Holes[i][j[0]], Holes[i][j[1]]]
                distances.append(distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path))
                for k in mutual_neighbors(edge, Matrix_out):
                    if k not in [vertex] + Matrix_out[vertex]:
                        outer_vertices = union(outer_vertices, k)
                    else:
                        for l in Holes[i]:
                            if k not in Matrix_out[l] and not is_line([k, l], other_edges) and not is_line([k, l],
                                                                                                         Holes[i]):
                                other_edges.append([k, l])
            for v in outer_vertices:
                distances.append(distance_vox(vertex, v, Centroids, coordinates, voxels, dist_path))
            for edge in other_edges:
                distances.append(distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path))
            Ind = distances.index(min(distances))
            if Ind < len(indices):
                edge = [Holes[i][item] for item in indices[Ind]]
            elif Ind < len(indices) + len(outer_vertices):
                edge = [vertex, outer_vertices[Ind - len(indices)]]
            else:
                edge = other_edges[Ind - len(indices) - len(outer_vertices)]
            if distances[Ind] > 3 * threshold:
                continue
            Register_Added.append(edge)
            Matrix_out = add_edge(edge[0], edge[1], Matrix_out)
            print("Added:", vertex, Holes[i], edge, distances[Ind])
            log_file_detailed.write("\nAdded: linkvertex {}, hole {}, edge {} with distance {}".format(
                vertex, Holes[i], edge, distances[Ind]))
            log_file_detailed.flush()
            out = False
            All = mutual_neighbors(edge, Matrix_out)
            for vertex_sub in All:
                computed_holes[vertex_sub], deg_3_edges_in_cross[vertex_sub], computed_open_edges[vertex_sub] = \
                    holes(vertex_sub, Matrix_out, list_boundary, Centroids, True, holes_degree=4)
            index -= 1
            break
    if out:
        break
print("")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

start_time = time()
print("\n\nCompute candidate edges to be added to close the holes...")
log_file_detailed.write("\n\nCompute candidate edges to be added to close the holes...")
log_file_detailed.flush()
log_file.write("\n\nCompute candidate edges to be added to close the holes...")
log_file.flush()
candidates_edges = [[] for item in range(num_pores)]
candidates_distances = [[] for item in range(num_pores)]
lst_compute = range(num_pores)
while lst_compute:
    vertex = lst_compute[0]
    candidates_edges[vertex] = []
    candidates_distances[vertex] = []
    print("\rCompute the candidates:", vertex, "out of", num_pores - 1, end="")
    if vertex in ignore:
        lst_compute = lst_compute[1:]
        continue
    Holes, edges_in_cross = [], []
    open_edges = open_edges_in_link(vertex, Matrix_out, Centroids)
    if len(open_edges) <= 75:
        Holes, edges_in_cross, _ = holes(vertex, Matrix_out, list_boundary, Centroids, True)
    elif vertex not in list_boundary:
        Holes = [unique(open_edges)]
    outer_vertices, other_edges = [], []
    for i in range(len(Holes)):
        priority = []
        for j in range(len(Holes[i])):
            if edges_in_cross and is_line([Holes[i][item] for item in [j, Next(j, 1, len(Holes[i]))]],
                                         edges_in_cross) and is_line(
                [Holes[i][item] for item in [j, Previous(j, 1, len(Holes[i]))]], edges_in_cross):
                priority.append(Holes[i][j])
        for j in range(len(Holes[i]) - 2):  # extract priority vertices
            if Holes[i][j] == Holes[i][j + 2]:
                priority = union(priority, Holes[i][j + 1])
        indices = []
        for j in combinations(range(len(Holes[i])), 2):
            j = list(j)
            edge = [Holes[i][item] for item in j]
            if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                continue
            if priority and not intersect(edge, priority) and len(Holes[i]) != 4:
                continue
            indices.append(j)
        for j in range(len(indices)):
            edge = [Holes[i][indices[j][0]], Holes[i][indices[j][1]]]
            if all([v in list_boundary_vox for v in Holes[i]]):
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                if dist > 2 * threshold:
                    continue
            if not is_line(edge, candidates_edges[vertex]):
                candidates_edges[vertex].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                candidates_distances[vertex].append(dist)
                for k in mutual_neighbors(edge, Matrix_out):
                    if k not in [vertex] + Matrix_out[vertex]:
                        outer_vertices = union(outer_vertices, k)
                    else:
                        for l in Holes[i]:
                            if k not in Matrix_out[l] and not is_line([k, l], other_edges) and not is_line([k, l],
                                                                                                         Holes[i]):
                                dist = distance_vox(k, l, Centroids, coordinates, voxels, dist_path)
                                if dist > 1.5 * threshold:
                                    continue
                                other_edges.append([k, l])
    for v in outer_vertices:
        edge = [vertex, v]
        if not is_line(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            candidates_distances[vertex].append(dist)
    for edge in other_edges:
        if not is_line(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            candidates_distances[vertex].append(dist)
    lst_compute = lst_compute[1:]
print("")

log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


start_time = time()
print("\n\nAdd short candidate edges to close the holes...")
log_file_detailed.write("\n\nAdd short candidate edges to close the holes...")
log_file_detailed.flush()
log_file.write("\n\nAdd short candidate edges to close the holes...")
log_file.flush()
count = 0
while True:
    minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(num_pores)]
    if min(minima) > 6 * threshold:
        break
    vertex = minima.index(min(minima))
    Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
    edge = candidates_edges[vertex][Ind]
    Register_Added.append(edge)
    Matrix_out = add_edge(edge[0], edge[1], Matrix_out)
    print("Added:", vertex, edge, candidates_distances[vertex][Ind])
    log_file_detailed.write(
        "\nAdded: linkvertex {}, edge {} with distance {}".format(vertex, edge, candidates_distances[vertex][Ind]))
    log_file_detailed.flush()
    All = union(edge, mutual_neighbors(edge, Matrix_out))
    for v_sub in All:
        candidates_edges[v_sub] = []
        candidates_distances[v_sub] = []
        Holes, edges_in_cross = [], []
        open_edges = open_edges_in_link(v_sub, Matrix_out, Centroids)
        if len(open_edges) <= 75:
            Holes, edges_in_cross, _ = holes(v_sub, Matrix_out, list_boundary, Centroids, True)
        elif v_sub not in list_boundary:
            Holes = [unique(open_edges)]
        outer_vertices, other_edges = [], []
        for i in range(len(Holes)):
            priority = []
            for j in range(len(Holes[i])):
                if edges_in_cross and is_line([Holes[i][item] for item in [j, Next(j, 1, len(Holes[i]))]],
                                             edges_in_cross) and is_line(
                    [Holes[i][item] for item in [j, Previous(j, 1, len(Holes[i]))]], edges_in_cross):
                    priority.append(Holes[i][j])
            for j in range(len(Holes[i]) - 2):  # extract priority vertices
                if Holes[i][j] == Holes[i][j + 2]:
                    priority = union(priority, Holes[i][j + 1])
            indices = []
            for j in combinations(range(len(Holes[i])), 2):
                j = list(j)
                edge = [Holes[i][item] for item in j]
                if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                    continue
                if priority and not intersect(edge, priority) and len(Holes[i]) != 4:
                    continue
                indices.append(j)
            for j in range(len(indices)):
                edge = [Holes[i][indices[j][0]], Holes[i][indices[j][1]]]
                if all([v in list_boundary_vox for v in Holes[i]]):
                    dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                    if dist > 2 * threshold:
                        continue
                if not is_line(edge, candidates_edges[v_sub]):
                    candidates_edges[v_sub].append(edge)
                    dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                    candidates_distances[v_sub].append(dist)
                    for k in mutual_neighbors(edge, Matrix_out):
                        if k not in [v_sub] + Matrix_out[v_sub]:
                            outer_vertices = union(outer_vertices, k)
                        else:
                            for l in Holes[i]:
                                if k not in Matrix_out[l] and not is_line([k, l], other_edges) and not is_line([k, l],
                                                                                                             Holes[i]):
                                    dist = distance_vox(k, l, Centroids, coordinates, voxels, dist_path)
                                    if dist > 1.5 * threshold:
                                        continue
                                    other_edges.append([k, l])
        for v in outer_vertices:
            edge = [v_sub, v]
            if not is_line(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                candidates_distances[v_sub].append(dist)
        for edge in other_edges:
            if not is_line(edge, candidates_edges[vertex]):
                candidates_edges[v_sub].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                candidates_distances[v_sub].append(dist)

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

start_time = time()
print("\n\nAdd edges to avoid football pores...")
log_file_detailed.write("\n\nAdd edges to avoid football pores...")
log_file_detailed.flush()
log_file.write("\n\nAdd edges to avoid football pores...")
log_file.flush()


# Add edges to avoid football pores. Start by collecting the candidates
print("Compute candidate edges to avoid football pores...")
K2 = []
for i in range(num_pores):
    Nei_i = Matrix_out[i]
    for j in [item for item in Nei_i if item > i]:
        K2.append([i, j])
candidates_edges = [[] for item in range(len(K2))]
candidates_distances = [[] for item in range(len(K2))]
lst_compute = range(len(K2))
count = 0
for index in lst_compute:
    candidate_edge = K2[index]
    candidates_edges[index] = []
    candidates_distances[index] = []
    count += 1
    print("\rComputing", count, "out of", len(K2), end="")
    if sublist(candidate_edge, list_boundary):
        continue
    lst = mutual_neighbors(candidate_edge, Matrix_out)
    if len(lst) > 6:
        continue
    if not cycles(edges_of(mutual_neighbors(candidate_edge, Matrix_out), Matrix_out), minimum_allowed=3):
        for w in set_diff(Matrix_out[candidate_edge[0]], [candidate_edge[1]] + Matrix_out[candidate_edge[1]]):
            candidates_edges[index].append([candidate_edge[1], w])
            dist = distance_vox(candidate_edge[1], w, Centroids, coordinates, voxels, dist_path)
            candidates_distances[index].append(dist)
        for w in set_diff(Matrix_out[candidate_edge[1]], [candidate_edge[0]] + Matrix_out[candidate_edge[0]]):
            candidates_edges[index].append([candidate_edge[0], w])
            dist = distance_vox(candidate_edge[0], w, Centroids, coordinates, voxels, dist_path)
            candidates_distances[index].append(dist)
        for edge in combinations(lst, 2):
            if edge[0] not in Matrix_out[edge[1]]:
                candidates_edges[index].append(list(edge))
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                candidates_distances[index].append(dist)


# Add short candidate edges...
print("\nAdd short candidate edges...")
while True:
    minima = [min(candidates_distances[index], default=math.inf) for index in range(len(K2))]
    if min(minima) > 6 * threshold:
        break
    index_edge = minima.index(min(minima))
    Ind = candidates_distances[index_edge].index(min(candidates_distances[index_edge]))
    edge = candidates_edges[index_edge][Ind]
    Register_Added.append(edge)
    Matrix_out = add_edge(edge[0], edge[1], Matrix_out)
    print("Added:", K2[index_edge], edge, candidates_distances[index_edge][Ind])
    log_file_detailed.write(
        "\nAdded to avoid a football pore between {} and {}: edge {} with distance {}".format(
            K2[index_edge][0], K2[index_edge][1], edge, candidates_distances[index_edge][Ind]))
    log_file_detailed.flush()
    K2.append(edge)
    candidates_edges.append([])
    candidates_distances.append([])
    All = union(edge, mutual_neighbors(edge, Matrix_out))
    All_ind = unique([ind(item, K2) for item in All])
    for ind_sub in All_ind:
        candidates_edges[ind_sub] = []
        candidates_distances[ind_sub] = []
        candidate_edge = K2[ind_sub]
        if sublist(candidate_edge, list_boundary):
            continue
        lst = mutual_neighbors(candidate_edge, Matrix_out)
        if len(lst) > 6:
            continue
        if not cycles(edges_of(mutual_neighbors(candidate_edge, Matrix_out), Matrix_out), minimum_allowed=3):
            for w in set_diff(Matrix_out[candidate_edge[0]], [candidate_edge[1]] + Matrix_out[candidate_edge[1]]):
                candidates_edges[ind_sub].append([candidate_edge[1], w])
                dist = distance_vox(candidate_edge[1], w, Centroids, coordinates, voxels, dist_path)
                candidates_distances[ind_sub].append(dist)
            for w in set_diff(Matrix_out[candidate_edge[1]], [candidate_edge[0]] + Matrix_out[candidate_edge[0]]):
                candidates_edges[ind_sub].append([candidate_edge[0], w])
                dist = distance_vox(candidate_edge[0], w, Centroids, coordinates, voxels, dist_path)
                candidates_distances[ind_sub].append(dist)
            for edge in combinations(lst, 2):
                if edge[0] not in Matrix_out[edge[1]]:
                    candidates_edges[ind_sub].append(list(edge))
                    dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                    candidates_distances[ind_sub].append(dist)

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


start_time = time()
print("\n\nAdd edges of degree 3...")
log_file_detailed.write("\n\nAdd edges of degree 3...")
log_file_detailed.flush()
log_file.write("\n\nAdd edges of degree 3...")
log_file.flush()
# Add edges of degree 3
for vertex in range(num_pores):
    link = vertexlink(vertex, Matrix_out)
    for tri in link:
        lst = set_diff(mutual_neighbors(tri, Matrix_out), [vertex] + Matrix_out[vertex])
        for outer_vertex in [item for item in lst if item > vertex]:
            dist = distance_vox(vertex, outer_vertex, Centroids, coordinates, voxels, dist_path)
            if dist < 3 * threshold:
                Matrix_out = add_edge(vertex, outer_vertex, Matrix_out)
                Register_Added.append([vertex, outer_vertex])
                print("Added around a triangle:", vertex, tri, [vertex, outer_vertex], dist)
                log_file_detailed.write(
                    "\nAdded around the triangle {}, edge {} with distance {}".format(tri, [vertex, outer_vertex],
                                                                                      dist))
                log_file_detailed.flush()

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


start_time = time()
print("\n\nCompute candidate edges within pseudoholes to be added...")
log_file_detailed.write("\n\nCompute candidate edges within pseudoholes to be added...")
log_file_detailed.flush()
log_file.write("\n\nCompute candidate edges within pseudoholes to be added...")
log_file.flush()
candidates_edges = [[] for item in range(num_pores)]
candidates_distances = [[] for item in range(num_pores)]
lst_compute = range(num_pores)
while lst_compute:
    vertex = lst_compute[0]
    candidates_edges[vertex] = []
    candidates_distances[vertex] = []
    print("\rCompute the candidates from pseudoholes:", vertex, "out of", num_pores - 1, end="")
    if vertex in ignore:
        lst_compute = lst_compute[1:]
        continue
    Holes = []
    open_edges = open_edges_in_link(vertex, Matrix_out, Centroids)
    if len(open_edges) <= 75:
        Holes = pseudoholes(vertex, Matrix_out, list_boundary, Centroids) + \
                holes(vertex, Matrix_out, list_boundary, Centroids, True)[0]
    elif vertex not in list_boundary:
        Holes = [unique(open_edges)]
    outer_vertices, other_edges = [], []
    for i in range(len(Holes)):
        indices = []
        for j in combinations(range(len(Holes[i])), 2):
            j = list(j)
            edge = [Holes[i][item] for item in j]
            if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                continue
            indices.append(j)
        for j in range(len(indices)):
            edge = [Holes[i][indices[j][0]], Holes[i][indices[j][1]]]
            if all([v in list_boundary_vox for v in Holes[i]]):
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                if dist > 2 * threshold:
                    continue
            if not is_line(edge, candidates_edges[vertex]):
                candidates_edges[vertex].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                candidates_distances[vertex].append(dist)
                for k in mutual_neighbors(edge, Matrix_out):
                    if k not in [vertex] + Matrix_out[vertex]:
                        outer_vertices = union(outer_vertices, k)
                    else:
                        for l in Holes[i]:
                            if k not in Matrix_out[l] and not is_line([k, l], other_edges) and not is_line([k, l],
                                                                                                         Holes[i]):
                                dist = distance_vox(k, l, Centroids, coordinates, voxels, dist_path)
                                if dist > 1.5 * threshold:
                                    continue
                                other_edges.append([k, l])
    for v in outer_vertices:
        edge = [vertex, v]
        if not is_line(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            candidates_distances[vertex].append(dist)
    for edge in other_edges:
        if not is_line(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            candidates_distances[vertex].append(dist)
    lst_compute = lst_compute[1:]
print("")


log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


start_time = time()
print("\n\nAdd short candidate edges to close the pseudoholes...")
log_file_detailed.write("\n\nAdd short candidate edges to close the pseudoholes...")
log_file_detailed.flush()
log_file.write("\n\nAdd short candidate edges to close the pseudoholes...")
log_file.flush()
count = 0
while True:
    minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(num_pores)]
    if min(minima) > 5 * threshold:
        break
    vertex = minima.index(min(minima))
    Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
    edge = candidates_edges[vertex][Ind]
    Register_Added.append(edge)
    Matrix_out = add_edge(edge[0], edge[1], Matrix_out)
    print("Added:", vertex, edge, candidates_distances[vertex][Ind])
    log_file_detailed.write(
        "\nAdded: linkvertex {}, edge {} with distance {}".format(vertex, edge, candidates_distances[vertex][Ind]))
    log_file_detailed.flush()
    All = union(edge, mutual_neighbors(edge, Matrix_out))
    for v_sub in All:
        candidates_edges[v_sub] = []
        candidates_distances[v_sub] = []
        Holes = []
        open_edges = open_edges_in_link(v_sub, Matrix_out, Centroids)
        if len(open_edges) <= 75:
            Holes = pseudoholes(v_sub, Matrix_out, list_boundary, Centroids) + \
                    holes(v_sub, Matrix_out, list_boundary, Centroids, True)[0]
        elif v_sub not in list_boundary:
            Holes = [unique(open_edges)]
        outer_vertices, other_edges = [], []
        for i in range(len(Holes)):
            indices = []
            for j in combinations(range(len(Holes[i])), 2):
                j = list(j)
                edge = [Holes[i][item] for item in j]
                if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                    continue
                indices.append(j)
            for j in range(len(indices)):
                edge = [Holes[i][indices[j][0]], Holes[i][indices[j][1]]]
                if all([v in list_boundary_vox for v in Holes[i]]):
                    dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                    if dist > 2 * threshold:
                        continue
                if not is_line(edge, candidates_edges[v_sub]):
                    candidates_edges[v_sub].append(edge)
                    dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                    candidates_distances[v_sub].append(dist)
                    for k in mutual_neighbors(edge, Matrix_out):
                        if k not in [v_sub] + Matrix_out[v_sub]:
                            outer_vertices = union(outer_vertices, k)
                        else:
                            for l in Holes[i]:
                                if k not in Matrix_out[l] and not is_line([k, l], other_edges) and not is_line([k, l],
                                                                                                             Holes[i]):
                                    dist = distance_vox(k, l, Centroids, coordinates, voxels, dist_path)
                                    if dist > 1.5 * threshold:
                                        continue
                                    other_edges.append([k, l])
        for v in outer_vertices:
            edge = [v_sub, v]
            if not is_line(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                candidates_distances[v_sub].append(dist)
        for edge in other_edges:
            if not is_line(edge, candidates_edges[vertex]):
                candidates_edges[v_sub].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                candidates_distances[v_sub].append(dist)


avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

log_file_detailed.write("\n\nRemove excessive edges...")
log_file_detailed.flush()
log_file.write("\nRemove excessive edges...")
log_file.flush()
print("\n\nRemove excessive edges...")

for vertex in range(num_pores):
    print("\rProgress:", vertex, "out of", num_pores - 1, end="")
    for triangle in vertexlink(vertex, Matrix_out):
        for edge_test in edges_vertexlink(vertex, Matrix_out):
            if intersect(edge_test, triangle):
                continue
            A, B = barycentric_intersection_in_line(triangle, edge_test[0], edge_test[1], Centroids)
            if min(A, B) >= 0.3:
                a, b, c = barycentric_intersection_in_triangle(triangle, edge_test, Centroids)
                if min(a, b, c) > 0:
                    dist = distance_vox(edge_test[0], edge_test[1], Centroids, coordinates, voxels, dist_path)
                    dists = [distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path) for edge in
                             combinations(triangle, 2)]
                    if dist > max([item + 0.5 * threshold for item in dists] + [2 * threshold]):
                        Register_Added = remove_index(Register_Added, ind(edge_test, Register_Added))
                        Register_Removed.append(edge_test)
                        Matrix_out = remove_edge(edge_test[0], edge_test[1], Matrix_out)
                        print("Removed:", vertex, triangle, edge_test, dist, dists)
                        log_file_detailed.write(
                            "\nRemoved excessive edge {} blocked by {} in the link of {}".format(edge_test, triangle,
                                                                                                 vertex))
                        log_file_detailed.flush()

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


# Collect candidates to close all remaining holes...
start_time = time()
print("\n\nCollect candidates to close all remaining holes...")
log_file_detailed.write("\n\nCollect candidates to close all remaining holes...")
log_file_detailed.flush()
log_file.write("\n\nCollect candidates to close all remaining holes...")
log_file.flush()
candidates_edges = [[] for item in range(num_pores)]
candidates_distances = [[] for item in range(num_pores)]
lst_compute = range(num_pores)
while lst_compute:
    vertex = lst_compute[0]
    candidates_edges[vertex] = []
    candidates_distances[vertex] = []
    print("\rCompute the candidates:", vertex, "out of", num_pores - 1, end="")
    if vertex in union(list_boundary, list_boundary_vox, ignore):
        lst_compute = lst_compute[1:]
        continue
    Holes = holes(vertex, Matrix_out, list_boundary, Centroids, True)[0]
    Matrix_reduced = reduced_vertexlink(vertex, Matrix_out, list_boundary, Centroids, return_Matrix=True)[1]
    Holes += [hole for hole in holes(vertex, Matrix_reduced, list_boundary, Centroids, True)[0] if
              not is_line_exact(hole, Holes)[0]]
    scale = 1
    if not Holes:
        Holes = pseudoholes(vertex, Matrix_out, list_boundary, Centroids)
        scale = 0.5
    outer_vertices, other_edges = [], []
    for i in range(len(Holes)):
        boundary_hole = len(intersect(Holes[i], union(list_boundary, list_boundary_vox))) >= len(Holes[i]) - 1
        for j in combinations(range(len(Holes[i])), 2):
            j = list(j)
            edge = [Holes[i][item] for item in j]
            if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                continue
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            if dist > scale * 10 * threshold:
                continue
            if all([v in list_boundary_vox for v in Holes[i]]):
                if dist > 3 * threshold:
                    continue
            if boundary_hole:
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                if dist > 6 * threshold:
                    continue
            if not is_line(edge, candidates_edges[vertex]):
                candidates_edges[vertex].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                candidates_distances[vertex].append(dist)
                for k in mutual_neighbors(edge, Matrix_out):
                    if k not in [vertex] + Matrix_out[vertex]:
                        outer_vertices = union(outer_vertices, k)
                    else:
                        for l in Holes[i]:
                            if k not in Matrix_out[l] and not is_line([k, l], other_edges) and not is_line([k, l],
                                                                                                         Holes[i]):
                                dist = distance_vox(k, l, Centroids, coordinates, voxels, dist_path)
                                if dist > 1.5 * threshold:
                                    continue
                                other_edges.append([k, l])
    for v in outer_vertices:
        edge = [vertex, v]
        dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
        if dist > scale * 10 * threshold:
            continue
        if not is_line(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            candidates_distances[vertex].append(dist)
    for edge in other_edges:
        dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
        if dist > scale * 10 * threshold:
            continue
        if not is_line(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            candidates_distances[vertex].append(dist)
    lst_compute = lst_compute[1:]
print("")

log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


start_time = time()
print("\n\nAdd short candidate edges to close remaining holes...")
log_file_detailed.write("\n\nAdd short candidate edges to close remaining holes...")
log_file_detailed.flush()
log_file.write("\n\nAdd short candidate edges to close remaining holes...")
log_file.flush()
while True:
    minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(num_pores)]
    if min(minima) > 10 * threshold:
        break
    vertex = minima.index(min(minima))
    Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
    edge = candidates_edges[vertex][Ind]
    Register_Added.append(edge)
    Matrix_out = add_edge(edge[0], edge[1], Matrix_out)
    print("Added:", vertex, edge, candidates_distances[vertex][Ind])
    log_file_detailed.write(
        "\nAdded: linkvertex {}, edge {} with distance {}".format(vertex, edge, candidates_distances[vertex][Ind]))
    log_file_detailed.flush()
    All = union(edge, mutual_neighbors(edge, Matrix_out))
    for v_sub in All:
        if v_sub in union(list_boundary, list_boundary_vox) or v_sub in ignore:
            continue
        candidates_edges[v_sub] = []
        candidates_distances[v_sub] = []
        Holes = holes(v_sub, Matrix_out, list_boundary, Centroids, True)[0]
        Matrix_temp = reduced_vertexlink(v_sub, Matrix_out, list_boundary, Centroids, return_Matrix=True)[1]
        Holes += [hole for hole in holes(v_sub, Matrix_temp, list_boundary, Centroids, True)[0] if
                  not is_line_exact(hole, Holes)[0]]
        scale = 1
        if not Holes:
            Holes = pseudoholes(v_sub, Matrix_out, list_boundary, Centroids)
            scale = 0.5
        outer_vertices, other_edges = [], []
        for i in range(len(Holes)):
            boundary_hole = len(intersect(Holes[i], union(list_boundary, list_boundary_vox))) >= len(Holes[i]) - 1
            for j in combinations(range(len(Holes[i])), 2):
                j = list(j)
                edge = [Holes[i][item] for item in j]
                if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                    continue
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                if dist > scale * 10 * threshold:
                    continue
                if all([v in list_boundary_vox for v in Holes[i]]):
                    if dist > 3 * threshold:
                        continue
                if boundary_hole:
                    dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                    if dist > 6 * threshold:
                        continue
                if not is_line(edge, candidates_edges[v_sub]):
                    candidates_edges[v_sub].append(edge)
                    dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                    candidates_distances[v_sub].append(dist)
                    for k in mutual_neighbors(edge, Matrix_out):
                        if k not in [v_sub] + Matrix_out[v_sub]:
                            outer_vertices = union(outer_vertices, k)
                        else:
                            for l in Holes[i]:
                                if k not in Matrix_out[l] and not is_line([k, l], other_edges) and not is_line([k, l],
                                                                                                             Holes[i]):
                                    dist = distance_vox(k, l, Centroids, coordinates, voxels, dist_path)
                                    if dist > 1.5 * threshold:
                                        continue
                                    other_edges.append([k, l])
        for v in outer_vertices:
            edge = [v_sub, v]
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            if dist > scale * 10 * threshold:
                continue
            if not is_line(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                candidates_distances[v_sub].append(dist)
        for edge in other_edges:
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            if dist > scale * 10 * threshold:
                continue
            if not is_line(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                candidates_distances[v_sub].append(dist)

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


# Collect candidates to close all remaining boundary holes...
start_time = time()
print("\n\nCollect candidates to close all remaining boundary holes...")
log_file_detailed.write("\n\nCollect candidates to close all remaining boundary holes...")
log_file_detailed.flush()
log_file.write("\n\nCollect candidates to close all remaining boundary holes...")
log_file.flush()
candidates_edges = [[] for item in range(num_pores)]
candidates_distances = [[] for item in range(num_pores)]
lst_compute = range(num_pores)
while lst_compute:
    vertex = lst_compute[0]
    candidates_edges[vertex] = []
    candidates_distances[vertex] = []
    print("\rCompute the candidates:", vertex, "out of", num_pores - 1, end="")
    if vertex in union(list_boundary, list_boundary_vox, ignore):
        lst_compute = lst_compute[1:]
        continue
    Holes = holes(vertex, Matrix_out, list_boundary, Centroids, True)[0]
    Matrix_reduced = reduced_vertexlink(vertex, Matrix_out, list_boundary, Centroids, return_Matrix=True)[1]
    Holes += [hole for hole in holes(vertex, Matrix_reduced, list_boundary, Centroids, True)[0] if
              not is_line_exact(hole, Holes)[0]]
    outer_vertices, other_edges = [], []
    for i in range(len(Holes)):
        for j in combinations(range(len(Holes[i])), 2):
            j = list(j)
            edge = [Holes[i][item] for item in j]
            if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                continue
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            if dist > 10 * threshold:
                continue
            if not is_line(edge, candidates_edges[vertex]):
                candidates_edges[vertex].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                candidates_distances[vertex].append(dist)
                for k in mutual_neighbors(edge, Matrix_out):
                    if k not in [vertex] + Matrix_out[vertex]:
                        outer_vertices = union(outer_vertices, k)
                    else:
                        for l in Holes[i]:
                            if k not in Matrix_out[l] and not is_line([k, l], other_edges) and not is_line([k, l],
                                                                                                         Holes[i]):
                                dist = distance_vox(k, l, Centroids, coordinates, voxels, dist_path)
                                if dist > 1.5 * threshold:
                                    continue
                                other_edges.append([k, l])
    for v in outer_vertices:
        edge = [vertex, v]
        dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
        if dist > 10 * threshold:
            continue
        if not is_line(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            candidates_distances[vertex].append(dist)
    for edge in other_edges:
        dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
        if dist > 10 * threshold:
            continue
        if not is_line(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            candidates_distances[vertex].append(dist)
    lst_compute = lst_compute[1:]
print("")

log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


start_time = time()
print("\n\nAdd short candidate edges to close remaining boundary holes...")
log_file_detailed.write("\n\nAdd short candidate edges to close remaining boundary holes...")
log_file_detailed.flush()
log_file.write("\n\nAdd short candidate edges to close remaining boundary holes...")
log_file.flush()
while True:
    minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(num_pores)]
    if min(minima) > 10 * threshold:
        break
    vertex = minima.index(min(minima))
    Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
    edge = candidates_edges[vertex][Ind]
    Register_Added.append(edge)
    Matrix_out = add_edge(edge[0], edge[1], Matrix_out)
    print("Added:", vertex, edge, candidates_distances[vertex][Ind])
    log_file_detailed.write(
        "\nAdded: linkvertex {}, edge {} with distance {}".format(vertex, edge, candidates_distances[vertex][Ind]))
    log_file_detailed.flush()
    All = union(edge, mutual_neighbors(edge, Matrix_out))
    for v_sub in All:
        if v_sub in union(list_boundary, list_boundary_vox) or v_sub in ignore:
            continue
        candidates_edges[v_sub] = []
        candidates_distances[v_sub] = []
        Holes = holes(v_sub, Matrix_out, list_boundary, Centroids, True)[0]
        Matrix_temp = reduced_vertexlink(v_sub, Matrix_out, list_boundary, Centroids, return_Matrix=True)[1]
        Holes += [hole for hole in holes(v_sub, Matrix_temp, list_boundary, Centroids, True)[0] if
                  not is_line_exact(hole, Holes)[0]]
        outer_vertices, other_edges = [], []
        for i in range(len(Holes)):
            for j in combinations(range(len(Holes[i])), 2):
                j = list(j)
                edge = [Holes[i][item] for item in j]
                if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                    continue
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                if dist > 10 * threshold:
                    continue
                if not is_line(edge, candidates_edges[v_sub]):
                    candidates_edges[v_sub].append(edge)
                    dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                    candidates_distances[v_sub].append(dist)
                    for k in mutual_neighbors(edge, Matrix_out):
                        if k not in [v_sub] + Matrix_out[v_sub]:
                            outer_vertices = union(outer_vertices, k)
                        else:
                            for l in Holes[i]:
                                if k not in Matrix_out[l] and not is_line([k, l], other_edges) and not is_line([k, l],
                                                                                                             Holes[i]):
                                    dist = distance_vox(k, l, Centroids, coordinates, voxels, dist_path)
                                    if dist > 1.5 * threshold:
                                        continue
                                    other_edges.append([k, l])
        for v in outer_vertices:
            edge = [v_sub, v]
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            if dist > 10 * threshold:
                continue
            if not is_line(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                candidates_distances[v_sub].append(dist)
        for edge in other_edges:
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            if dist > 10 * threshold:
                continue
            if not is_line(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                candidates_distances[v_sub].append(dist)

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


# Detect football pores and collect candidates to repair non-Hamiltonian links...
start_time = time()
print("\n\nDetect football pores and collect candidates to repair non-Hamiltonian links...")
log_file_detailed.write("\n\nDetect football pores and collect candidates to repair non-Hamiltonian links...")
log_file_detailed.flush()
log_file.write("\n\nDetect football pores and collect candidates to repair non-Hamiltonian links...")
log_file.flush()
candidates_edges = [[] for item in range(num_pores)]
candidates_distances = [[] for item in range(num_pores)]
lst_Hamiltonian, lst_not_Hamiltonian = [], []
football_pores, football_pores_neighbors = [], []
loop_over = range(num_pores)
avoid = union(list_boundary, list_boundary_vox, ignore)
while loop_over:
    vertex = loop_over[0]
    print("\rCompute the candidates:", vertex, "out of", num_pores - 1, end="")
    if vertex in avoid:
        loop_over = loop_over[1:]
        continue
    Matrix_reduced, triangles = reduced_vertexlink(vertex, Matrix_out, list_boundary, Centroids, return_Matrix=True,
                                                   return_triangles=True)[1:]
    tets = [[vertex] + item for item in triangles]
    if is_Hamiltonian(vertex, Matrix_reduced, list_boundary, Centroids, predetermined_tetrahedra=tets)[1]:
        lst_Hamiltonian.append(vertex)
        loop_over = loop_over[1:]
        continue
    else:
        if is_Hamiltonian(vertex, Matrix_out, list_boundary, Centroids)[1]:
            lst_Hamiltonian.append(vertex)
            loop_over = loop_over[1:]
            continue
        else:
            lst_not_Hamiltonian.append(vertex)
    if len(Matrix_out[vertex]) < 6:
        loop_over = loop_over[1:]
        continue  # this vertex might represent a football pore
    if len(Matrix_out[vertex]) > 40:
        loop_over = loop_over[1:]
        continue
    lst = [len(mutual_neighbors([vertex, v], Matrix_out)) for v in Matrix_out[vertex]]
    Ind = sorted(range(len(lst)), key=lambda k: lst[k])
    # avoid repairing football pores
    neighbors_ordered = [Matrix_out[vertex][item] for item in Ind]
    digon_involved = False
    for v in neighbors_ordered:
        if len(Matrix_out[v]) > 5 or v in list_boundary:
            continue
        if cycles(edges_of(mutual_neighbors([vertex, v], Matrix_out), Matrix_out), minimum_allowed=3):
            continue
        Matrix_temp = remove_edge(v, vertex, Matrix_out)
        if is_Hamiltonian(vertex, Matrix_temp, list_boundary, Centroids)[1]:
            digon_involved = True
            football_pores.append(v)
            ignore.append(v)
            avoid.append(v)
            football_pore_neighbors = Matrix_out[v]
            football_pores_neighbors.append(football_pore_neighbors)
            for u in football_pore_neighbors:
                Matrix_out = remove_edge(u, v, Matrix_out)
            print("\nfootball pore found: {} - neighbors: {}".format(v, football_pore_neighbors))
            log_file_detailed.write("\nfootball pore found: {} - neighbors: {}".format(v, football_pore_neighbors))
            log_file_detailed.flush()
            for v_sub in football_pore_neighbors:
                if v_sub in union(list_boundary, list_boundary_vox, ignore):
                    continue
                Matrix_reduced, triangles = reduced_vertexlink(v_sub, Matrix_out, list_boundary, Centroids,
                                                               return_Matrix=True, return_triangles=True)[1:]
                lst_Hamiltonian = set_diff(lst_Hamiltonian, v_sub)
                lst_not_Hamiltonian = set_diff(lst_not_Hamiltonian, v_sub)
                loop_over = set_diff(loop_over, v_sub)
                tets = [[v_sub] + item for item in triangles]
                if is_Hamiltonian(v_sub, Matrix_reduced, list_boundary, Centroids, predetermined_tetrahedra=tets)[1]:
                    lst_Hamiltonian.append(v_sub)
                else:
                    if is_Hamiltonian(v_sub, Matrix_out, list_boundary, Centroids)[1]:
                        lst_Hamiltonian.append(v_sub)
                    else:
                        lst_not_Hamiltonian.append(v_sub)
                        loop_over.append(v_sub)
            break
    if digon_involved:
        continue
    second_neighbors = union(*[Matrix_out[item] for item in Matrix_out[vertex]])
    initial_Hamiltonian = intersect(second_neighbors, lst_Hamiltonian)
    candidates_Hamiltonian, candidates_add, candidates_add_dist = [], [], []
    for i in Matrix_out[vertex]:
        for j in [item for item in set_diff(Matrix_out[vertex], Matrix_out[i]) if item > i]:
            edge = [i, j]
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            if dist > 10 * threshold:
                continue
            Matrix_temp = add_edge(edge[0], edge[1], Matrix_out)
            Matrix_reduced, triangles = reduced_vertexlink(vertex, Matrix_temp, list_boundary, Centroids,
                                                           return_Matrix=True, return_triangles=True)[1:]
            tets = [[vertex] + item for item in triangles]
            check = is_Hamiltonian(vertex, Matrix_reduced, list_boundary, Centroids, predetermined_tetrahedra=tets)[1]
            if not check:
                check = is_Hamiltonian(vertex, Matrix_temp, list_boundary, Centroids)[1]
            if check:
                candidates_add.append(edge)
                candidates_add_dist.append(dist)
                candidates_Hamiltonian.append([vertex])
                for v in set_diff(second_neighbors, vertex):
                    if v in union(list_boundary, list_boundary_vox, ignore):
                        continue
                    if v in lst_Hamiltonian and v not in edge:
                        candidates_Hamiltonian[-1].append(v)
                        continue
                    if v not in union(edge, mutual_neighbors(edge, Matrix_out)):
                        continue
                    Matrix_reduced, triangles = reduced_vertexlink(v, Matrix_temp, list_boundary, Centroids,
                                                                   return_Matrix=True, return_triangles=True)[1:]
                    tets = [[v] + item for item in triangles]
                    if is_Hamiltonian(v, Matrix_reduced, list_boundary, Centroids, predetermined_tetrahedra=tets)[1]:
                        candidates_Hamiltonian[-1].append(v)
                        continue
                    if is_Hamiltonian(v, Matrix_temp, list_boundary, Centroids)[1]:
                        candidates_Hamiltonian[-1].append(v)
    diff = [len(item) - len(initial_Hamiltonian) for item in candidates_Hamiltonian]
    m = max(diff, default=0)
    if m > 0:
        Ind = [index for index in range(len(candidates_Hamiltonian)) if diff[index] == m]
        for index in set_diff(range(len(candidates_add_dist)), Ind):
            candidates_add_dist[index] = math.inf
        Ind = candidates_add_dist.index(min(candidates_add_dist))  # add the shortest edge
        edge = candidates_add[Ind]
        dist = candidates_add_dist[Ind]
        candidates_edges[vertex].append(edge)
        candidates_distances[vertex].append(dist)
    loop_over = loop_over[1:]
print("")


start_time = time()
print("\n\nAdd short candidate edges to repair non-Hamiltonian links...")
log_file_detailed.write("\n\nAdd short candidate edges to repair non-Hamiltonian links...")
log_file_detailed.flush()
log_file.write("\n\nAdd short candidate edges to repair non-Hamiltonian links...")
log_file.flush()
while True:
    minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(num_pores)]
    if min(minima) > 10 * threshold:
        break
    vertex = minima.index(min(minima))
    Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
    edge = candidates_edges[vertex][Ind]
    Register_Added.append(edge)
    Matrix_out = add_edge(edge[0], edge[1], Matrix_out)
    Register_Removed = remove_index(Register_Removed, ind(edge, Register_Removed))
    Register_Added.append(edge)
    dist = candidates_distances[vertex][Ind]
    print("Repair the link of {} by adding edge {} with length {}".format(vertex, edge, dist))
    log_file_detailed.write("\nRepair the link of {} by adding edge {} with length {}".format(vertex, edge, dist))
    log_file_detailed.flush()
    for v_sub in edge + mutual_neighbors(edge, Matrix_out):
        if v_sub in union(list_boundary, list_boundary_vox, ignore):
            continue
        candidates_edges[v_sub] = []
        candidates_distances[v_sub] = []
        lst_Hamiltonian = set_diff(lst_Hamiltonian, v_sub)
        lst_not_Hamiltonian = set_diff(lst_not_Hamiltonian, v_sub)
        Matrix_reduced, triangles = reduced_vertexlink(v_sub, Matrix_out, list_boundary, Centroids, return_Matrix=True,
                                                       return_triangles=True)[1:]
        tets = [[v_sub] + item for item in triangles]
        if is_Hamiltonian(v_sub, Matrix_reduced, list_boundary, Centroids, predetermined_tetrahedra=tets)[1]:
            lst_Hamiltonian = union(lst_Hamiltonian, v_sub)
            continue
        else:
            if is_Hamiltonian(v_sub, Matrix_out, list_boundary, Centroids)[1]:
                lst_Hamiltonian = union(lst_Hamiltonian, v_sub)
                continue
            else:
                lst_not_Hamiltonian = union(lst_not_Hamiltonian, v_sub)
        if len(Matrix_out[v_sub]) < 6:
            loop_over = loop_over[1:]
            continue  # this vertex might represent a football pore
        if len(Matrix_out[v_sub]) > 40:
            loop_over = loop_over[1:]
            continue
        second_neighbors = union(*[Matrix_out[item] for item in Matrix_out[v_sub]])
        initial_Hamiltonian = intersect(second_neighbors, lst_Hamiltonian)
        candidates_Hamiltonian, candidates_add, candidates_add_dist = [], [], []
        for i in Matrix_out[v_sub]:
            for j in [item for item in set_diff(Matrix_out[v_sub], Matrix_out[i]) if item > i]:
                edge = [i, j]
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                if dist > 10 * threshold:
                    continue
                Matrix_temp = add_edge(edge[0], edge[1], Matrix_out)
                Matrix_reduced, triangles = reduced_vertexlink(v_sub, Matrix_temp, list_boundary, Centroids,
                                                               return_Matrix=True, return_triangles=True)[1:]
                tets = [[v_sub] + item for item in triangles]
                check = is_Hamiltonian(v_sub, Matrix_reduced, list_boundary, Centroids,
                                       predetermined_tetrahedra=tets)[1]
                if not check:
                    check = is_Hamiltonian(v_sub, Matrix_temp, list_boundary, Centroids)[1]
                if check:
                    candidates_add.append(edge)
                    candidates_add_dist.append(dist)
                    candidates_Hamiltonian.append([v_sub])
                    for v in set_diff(second_neighbors, v_sub):
                        if v in union(list_boundary, list_boundary_vox, ignore):
                            continue
                        if v in lst_Hamiltonian and v not in edge:
                            candidates_Hamiltonian[-1].append(v)
                            continue
                        if v not in union(edge, mutual_neighbors(edge, Matrix_out)):
                            continue
                        Matrix_reduced, triangles = reduced_vertexlink(v, Matrix_temp, list_boundary, Centroids,
                                                                       return_Matrix=True, return_triangles=True)[1:]
                        tets = [[v] + item for item in triangles]
                        if is_Hamiltonian(
                                v, Matrix_reduced, list_boundary, Centroids, predetermined_tetrahedra=tets)[1]:
                            candidates_Hamiltonian[-1].append(v)
                            continue
                        if is_Hamiltonian(v, Matrix_temp, list_boundary, Centroids)[1]:
                            candidates_Hamiltonian[-1].append(v)
        diff = [len(item) - len(initial_Hamiltonian) for item in candidates_Hamiltonian]
        m = max(diff, default=0)
        if m > 0:
            Ind = [index for index in range(len(candidates_Hamiltonian)) if diff[index] == m]
            for index in set_diff(range(len(candidates_add_dist)), Ind):
                candidates_add_dist[index] = math.inf
            Ind = candidates_add_dist.index(min(candidates_add_dist))  # add the shortest edge
            edge = candidates_add[Ind]
            dist = candidates_add_dist[Ind]
            candidates_edges[v_sub].append(edge)
            candidates_distances[v_sub].append(dist)

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
avg_deg += sum([len(item) for item in football_pores_neighbors])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


# Add exceptionally long edges...
start_time = time()
print("\n\nAdd exceptionally long edges...")
log_file_detailed.write("\n\nAdd exceptionally long edges...")
log_file_detailed.flush()
log_file.write("\n\nAdd exceptionally long edges...")
log_file.flush()
print("Compute the candidate edges...")
log_file_detailed.write("\nCompute the candidate edges...")
log_file_detailed.flush()
log_file.write("\nCompute the candidate edges...")
log_file.flush()

candidates_edges = [[] for item in range(num_pores)]
candidates_distances = [[] for item in range(num_pores)]
lst_compute = range(num_pores)
while lst_compute:
    vertex = lst_compute[0]
    candidates_edges[vertex] = []
    candidates_distances[vertex] = []
    print("\rCompute the candidates:", vertex, "out of", num_pores - 1, end="")
    if vertex in union(list_boundary, list_boundary_vox, ignore):
        lst_compute = lst_compute[1:]
        continue
    Holes = holes(vertex, Matrix_out, list_boundary, Centroids, True)[0]
    Matrix_reduced = reduced_vertexlink(vertex, Matrix_out, list_boundary, Centroids, return_Matrix=True)[1]
    Holes += [hole for hole in holes(vertex, Matrix_reduced, list_boundary, Centroids, True)[0] if
              not is_line_exact(hole, Holes)[0]]
    scale = 1
    if not Holes:
        Holes = pseudoholes(vertex, Matrix_out, list_boundary, Centroids)
        scale = 0.25
    outer_vertices, other_edges = [], []
    for i in range(len(Holes)):
        for j in combinations(range(len(Holes[i])), 2):
            j = list(j)
            edge = [Holes[i][item] for item in j]
            if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                continue
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            if dist > scale * 20 * threshold:
                continue
            if not is_line(edge, candidates_edges[vertex]):
                candidates_edges[vertex].append(edge)
                candidates_distances[vertex].append(dist)
                for k in mutual_neighbors(edge, Matrix_out):
                    if k not in [vertex] + Matrix_out[vertex]:
                        outer_vertices = union(outer_vertices, k)
                    else:
                        for l in Holes[i]:
                            if k not in Matrix_out[l] and not is_line([k, l], other_edges) and not is_line(
                                    [k, l], Holes[i]):
                                dist = distance_vox(k, l, Centroids, coordinates, voxels, dist_path)
                                if dist > 1.5 * threshold:
                                    continue
                                other_edges.append([k, l])
    for v in outer_vertices:
        edge = [vertex, v]
        dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
        if dist > scale * 20 * threshold:
            continue
        if not is_line(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            candidates_distances[vertex].append(dist)
    for edge in other_edges:
        dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
        if dist > scale * 20 * threshold:
            continue
        if not is_line(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            candidates_distances[vertex].append(dist)
    lst_compute = lst_compute[1:]
print("")

start_time = time()
print("Add (probably long) candidate edges to close remaining holes...")
log_file_detailed.write("\nAdd (probably long) candidate edges to close remaining holes...")
log_file_detailed.flush()
log_file.write("\nAdd (probably long) candidate edges to close remaining holes...")
log_file.flush()
count = 0
while True:
    minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(num_pores)]
    if min(minima) > 20 * threshold:
        break
    vertex = minima.index(min(minima))
    Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
    edge = candidates_edges[vertex][Ind]
    Register_Removed = remove_index(Register_Removed, ind(edge, Register_Removed))
    Register_Added.append(edge)
    Matrix_out = add_edge(edge[0], edge[1], Matrix_out)
    print("Added:", vertex, edge, candidates_distances[vertex][Ind])
    log_file_detailed.write("\nAdded: linkvertex {}, edge {} with distance {}".format(
        vertex, edge, candidates_distances[vertex][Ind]))
    log_file_detailed.flush()
    All = union(edge, mutual_neighbors(edge, Matrix_out))
    for v_sub in All:
        if v_sub in union(list_boundary, list_boundary_vox) or vertex in ignore:
            continue
        candidates_edges[v_sub] = []
        candidates_distances[v_sub] = []
        Holes = holes(v_sub, Matrix_out, list_boundary, Centroids, True)[0]
        Matrix_reduced = reduced_vertexlink(v_sub, Matrix_out, list_boundary, Centroids, return_Matrix=True)[1]
        Holes += [hole for hole in holes(v_sub, Matrix_reduced, list_boundary, Centroids, True)[0] if
                  not is_line_exact(hole, Holes)[0]]
        scale = 1
        if not Holes:
            Holes = pseudoholes(v_sub, Matrix_out, list_boundary, Centroids)
            scale = 0.25
        outer_vertices, other_edges = [], []
        for i in range(len(Holes)):
            for j in combinations(range(len(Holes[i])), 2):
                j = list(j)
                edge = [Holes[i][item] for item in j]
                if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                    continue
                dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
                if dist > scale * 20 * threshold:
                    continue
                if not is_line(edge, candidates_edges[v_sub]):
                    candidates_edges[v_sub].append(edge)
                    candidates_distances[v_sub].append(dist)
                    for k in mutual_neighbors(edge, Matrix_out):
                        if k not in [v_sub] + Matrix_out[v_sub]:
                            outer_vertices = union(outer_vertices, k)
                        else:
                            for l in Holes[i]:
                                if k not in Matrix_out[l] and not is_line([k, l], other_edges) and not is_line([k, l],
                                                                                                             Holes[i]):
                                    dist = distance_vox(k, l, Centroids, coordinates, voxels, dist_path)
                                    if dist > 1.5 * threshold:
                                        continue
                                    other_edges.append([k, l])
        for v in outer_vertices:
            edge = [v_sub, v]
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            if dist > scale * 20 * threshold:
                continue
            if not is_line(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                candidates_distances[v_sub].append(dist)
        for edge in other_edges:
            dist = distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path)
            if dist > scale * 20 * threshold:
                continue
            if not is_line(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                candidates_distances[v_sub].append(dist)

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
avg_deg += sum([len(item) for item in football_pores_neighbors])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

remaining_interior_holes = []
for vertex in range(num_pores):
    if vertex in union(list_boundary, list_boundary_vox, ignore):
        continue
    Holes = holes(vertex, Matrix_out, list_boundary, Centroids, True)[0]
    if Holes:
        remaining_interior_holes.append(vertex)
        continue
    Holes = pseudoholes(vertex, Matrix_out, list_boundary, Centroids)
    if Holes:
        remaining_interior_holes.append(vertex)
        continue
    Matrix_reduced = reduced_vertexlink(vertex, Matrix_out, list_boundary, Centroids, return_Matrix=True)[1]
    Holes = holes(vertex, Matrix_reduced, list_boundary, Centroids, True)[0]
    if Holes:
        remaining_interior_holes.append(vertex)
print("{} remaining holes".format(len(remaining_interior_holes)))
log_file_detailed.write("\n\n{} interior vertices might still have holes in their links".
                        format(len(remaining_interior_holes)))
log_file_detailed.flush()
log_file.write("\n{} interior vertices might still have holes in their links".format(len(remaining_interior_holes)))
log_file.flush()

log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

K6 = []
for i in range(num_pores):
    Nei_i = Matrix_out[i]
    for j in [item for item in Nei_i if item > i]:
        Nei_i_j = intersect(Nei_i, Matrix_out[j])
        for k in [item for item in Nei_i_j if item > j]:
            Nei_i_j_k = intersect(Nei_i_j, Matrix_out[k])
            for l in [item for item in Nei_i_j_k if item > k]:
                Nei_i_j_k_l = intersect(Nei_i_j_k, Matrix_out[l])
                for m in [item for item in Nei_i_j_k_l if item > l]:
                    Nei_i_j_k_l_m = intersect(Nei_i_j_k_l, Matrix_out[m])
                    for n in [item for item in Nei_i_j_k_l_m if item > m]:
                        K6.append([i, j, k, l, m, n])

# Remove excessive edges from 6-cliques...
start_time = time()
print("\n\nRemove excessive edges from 6-cliques...")
log_file_detailed.write("\n\nRemove excessive edges from 6-cliques...")
log_file_detailed.flush()
log_file.write("\n\nRemove excessive edges from 6-cliques...")
log_file.flush()
print("Compute the candidate edges...")
log_file_detailed.write("\nCompute the candidate edges...")
log_file_detailed.flush()
log_file.write("\nCompute the candidate edges...")
log_file.flush()

candidates_edges = [[] for item in K6]
candidates_distances = [[] for item in K6]
for i in range(len(K6)):
    print("\rCompute the candidates:", i + 1, "out of", len(K6), end="")
    lst = K6[i]
    Holes_all = []
    for v in lst:
        Holes = holes(v, Matrix_out, list_boundary, Centroids, True)[0]
        Matrix_reduced = reduced_vertexlink(v, Matrix_out, list_boundary, Centroids, return_Matrix=True)[1]
        Holes += [hole for hole in holes(v, Matrix_reduced, list_boundary,
                                         Centroids, True)[0] if not is_line_exact(hole, Holes)[0]]
        Holes_all.append(Holes)
    edges, distances = [], []
    for edge in combinations(lst, 2):
        lst_sub = mutual_neighbors(edge, Matrix_out)
        lst_sub = set_diff(lst_sub, union(list_boundary, list_boundary_vox))
        Matrix_temp = remove_edge(edge[0], edge[1], Matrix_out)
        removable = True
        for v in lst_sub:
            if is_Hamiltonian(v, Matrix_out, list_boundary, Centroids)[1] and not is_Hamiltonian(
                    v, Matrix_temp, list_boundary, Centroids)[1]:
                removable = False
                break
            for j in range(len(lst)):
                v = lst[j]
                if is_line(v, union(list_boundary, list_boundary_vox, ignore)):
                    continue
                Holes = holes(v, Matrix_temp, list_boundary, Centroids, True)[0]
                Matrix_reduced = reduced_vertexlink(v, Matrix_temp, list_boundary, Centroids, return_Matrix=True)[1]
                Holes += [hole for hole in holes(v, Matrix_reduced, list_boundary, Centroids, True)[0] if
                          not is_line_exact(hole, Holes)[0]]
                if any([not is_line_exact(hole, Holes_all[j])[0] for hole in Holes]):
                    removable = False
                    break
            if not removable:
                break
        if removable:
            edges.append(edge)
            distances.append(distance_vox(edge[0], edge[1], Centroids, coordinates, voxels, dist_path))
    if distances:
        Ind = distances.index(max(distances))
        if distances[Ind] < threshold:
            candidates_distances[i].append(-math.inf)
            continue
        candidates_edges[i].append(list(edges[Ind]))
        candidates_distances[i].append(distances[Ind])
    else:
        candidates_distances[i].append(-math.inf)
print("")


print("Remove excessive edges...")
log_file_detailed.write("\nRemove excessive edges...")
log_file_detailed.flush()
log_file.write("\nRemove excessive edges...")
log_file.flush()
count = 0
while True:
    lst = [candidates_distances[item][0] for item in range(len(candidates_distances))]
    m = max(lst)
    if m == -math.inf:
        break
    Ind = lst.index(m)
    edge = candidates_edges[Ind][0]
    Register_Added = remove_index(Register_Added, ind(edge, Register_Added))
    Register_Removed.append(edge)
    Matrix_out = remove_edge(edge[0], edge[1], Matrix_out)
    print("Removed from 6-clique:", K6[Ind], edge, candidates_distances[Ind][0])
    log_file_detailed.write("\nRemoved from 6-clique {} - edge {} with distance {}".format(
        K6[Ind], edge, candidates_distances[Ind][0]))
    log_file_detailed.flush()
    Ind = ind(edge, K6)
    for i in Ind:
        candidates_edges[i] = []
        candidates_distances[i] = [-math.inf]

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
avg_deg += sum([len(item) for item in football_pores_neighbors])
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


start_time = time()
print("\n\nCompute the list of interior vertices with repairable or proper links...")
log_file_detailed.write("\n\nCompute the list of interior vertices with repairable or proper links...")
log_file_detailed.flush()
log_file.write("\n\nCompute the list of interior vertices with repairable or proper links...")
log_file.flush()
lst_repairable = [False for item in range(num_pores)]
lst_proper = [False for item in range(num_pores)]
for vertex in range(num_pores):
    print("\rCompute the list of vertices with repairable or proper links:", vertex, "out of", num_pores - 1,
          end="")
    if repairable_link(vertex, Matrix_out, list_boundary, Centroids):
        lst_repairable[vertex] = True
    if is_vertexlink_proper(vertex, Matrix_out, list_boundary, Centroids):
        lst_proper[vertex] = True
print("")
print("# repairable:", len([item for item in list_interior if lst_repairable[item]]), "out of", len(list_interior))
print("# proper:", len([item for item in list_interior if lst_proper[item]]), "out of", len(list_interior))
log_file_detailed.write("\n# repairable: {} out of {}".format(len([item for item in list_interior
                                                                   if lst_repairable[item]]), len(list_interior)))
log_file_detailed.write("\n# proper: {} out of {}".format(len([item for item in list_interior
                                                               if lst_proper[item]]), len(list_interior)))
log_file_detailed.flush()
log_file.write("\n# repairable: {} out of {}".format(len([item for item in list_interior
                                                          if lst_repairable[item]]), len(list_interior)))
log_file.write("\n# proper: {} out of {}".format(len([item for item in list_interior if lst_proper[item]]),
                                                 len(list_interior)))
log_file.flush()


K5 = find_K5s(Matrix_out)
remaining = []
for vertex in [item for item in list_interior if not lst_repairable[item] and not lst_proper[item]]:
    if vertex not in union(list_boundary_vox, list_boundary):
        remaining.append(vertex)
        continue
    _, _, problems = is_triangulatable(vertex, Matrix_out, list_boundary, Centroids)
    if set_diff(problems, union(list_boundary_vox, list_boundary)):
        remaining.append(vertex)

log_file_detailed.write("\n\nChange tunnels of 5-cliques...")
log_file_detailed.flush()
log_file.write("\n\nChange tunnels of 5-cliques...")
log_file.flush()
print("changing tunnels...")
changed_tunnels = [[], []]
for vertex in remaining:
    K5_temp = [K5[item] for item in ind(vertex, K5)]
    T, _, problems = is_triangulatable(vertex, Matrix_out, list_boundary, Centroids, changed_tunnels)
    if T:
        continue
    K5_temp = [item for item in K5_temp if intersect(item, problems)]
    for K5_current in K5_temp:
        tunnels = [list(item) for item in combinations(K5_current, 2)]
        K5_current_sub = set_diff(K5_current, index_inf)
        diff = []
        tunnel_initial = tunnel_K5_geometric(K5_current, Centroids)
        Ind = ind(K5_current, changed_tunnels[0])
        if Ind:
            tunnel_initial = changed_tunnels[1][Ind[0]]
        T_repaired = is_triangulatable(K5_current_sub, Matrix_out, list_boundary, Centroids, changed_tunnels,
                                       try_tunnels=K5_current)[0]
        if tunnel_initial in tunnels:
            Ind = tunnels.index(tunnel_initial)
            initial_repaired = sum([T_repaired[item][Ind] for item in range(len(K5_current_sub))])
        else:
            initial_repaired = sum(is_triangulatable(K5_current_sub, Matrix_out, list_boundary, Centroids,
                                                     changed_tunnels)[0])
        T_repaired = np.array(T_repaired)
        diff = np.sum(T_repaired, axis=0) - initial_repaired
        m = max(diff)
        print(vertex, list(diff), m)
        if m <= 0:
            continue
        print(("Tunnel changed: vertex {}, 5-clique: {}, gains: {}  |  {}".format(vertex, K5_current, list(diff), m)))
        log_file_detailed.write("\nTunnel changed: vertex {}, 5-clique: {}, gains: {}  |  {}".format(vertex, K5_current,
                                                                                                     list(diff), m))
        log_file_detailed.flush()
        Ind_repaired = [item for item in range(len(diff)) if diff[item] == m]
        if len(Ind_repaired) == 1:
            tunnel = tunnels[Ind_repaired[0]]
        else:
            T, _, problems = is_triangulatable(K5_current_sub, Matrix_out, list_boundary, Centroids, changed_tunnels,
                                               try_tunnels=K5_current, with_repair=False)
            num_of_problems = [len(item) for item in problems]
            m = min(num_of_problems)
            Ind_not_repaired = [item for item in range(len(num_of_problems)) if num_of_problems[item] == m]
            Ind = intersect(Ind_repaired, Ind_not_repaired)
            if Ind:
                Ind = Ind[0]
            else:
                Ind = Ind_repaired[0]
            tunnel = tunnels[Ind]

        Ind = ind(K5_current, changed_tunnels[0])
        if Ind:
            Ind = Ind[0]
            changed_tunnels[1][Ind] = tunnel
        else:
            changed_tunnels[0].append(K5_current)
            changed_tunnels[1].append(tunnel)
        if is_triangulatable(vertex, Matrix_out, list_boundary, Centroids, changed_tunnels)[0]:
            break


log_file_detailed.write("\n\nCompute the list of interior vertices with repairable or proper links...")
log_file_detailed.flush()
log_file.write("\n\nCompute the list of interior vertices with repairable or proper links...")
log_file.flush()
lst_proper = [False for item in range(num_pores)]
for vertex in range(num_pores):
    print("\rCompute the list of vertices with proper links:", vertex, "out of", num_pores - 1, end="")
    if is_vertexlink_proper(vertex, Matrix_out, list_boundary, Centroids, changed_tunnels):
        lst_proper[vertex] = True
print("")
print("# repairable:", len([item for item in list_interior if lst_repairable[item]]), "out of", len(list_interior))
print("# proper:", len([item for item in list_interior if lst_proper[item]]), "out of", len(list_interior))
log_file_detailed.write("\n# repairable: {} out of {}".format(len([item for item in list_interior
                                                                   if lst_repairable[item]]), len(list_interior)))
log_file_detailed.write("\n# proper: {} out of {}".format(len([item for item in list_interior
                                                               if lst_proper[item]]), len(list_interior)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n# repairable: {} out of {}".format(len([item for item in list_interior
                                                          if lst_repairable[item]]), len(list_interior)))
log_file.write("\n# proper: {} out of {}".format(len([item for item in list_interior
                                                      if lst_proper[item]]), len(list_interior)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


###############################
# Reconstruct the triangulation
###############################
start_time = time()
print("Reconstruct the triangulation at a first try...")
log_file_detailed.write("\n\nReconstruct the triangulation at a first try...")
log_file_detailed.flush()
log_file.write("\n\nReconstruct the triangulation at a first try...")
log_file.flush()

triangulation_reconstructed = []
lst_repeat = []
for vertex in ordered_pores:
    initial_len = len(triangulation_reconstructed)
    T, triangulation_reconstructed, _ = is_triangulatable(vertex, Matrix_out, list_boundary, Centroids, changed_tunnels,
                                                          with_repair=False, predetermined=triangulation_reconstructed)
    if not T:
        triangulation_reconstructed = triangulation_reconstructed[:initial_len]
        lst_repeat.append(vertex)
for vertex in lst_repeat:
    _, triangulation_reconstructed, _ = is_triangulatable(vertex, Matrix_out, list_boundary, Centroids, changed_tunnels,
                                                          predetermined=triangulation_reconstructed)

edges_defect_link, edges_defect_link_interior = [], []
count, count_interior = 0, 0
for i in range(num_pores):
    Nei = Matrix_out[i]
    for j in [item for item in Nei if item > i]:
        if sublist([i, j], union(list_boundary, list_boundary_vox)):
            continue
        count += 1
        if not intersect([i, j], union(list_boundary, list_boundary_vox)):
            count_interior += 1
            if not is_edge_proper_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])
                edges_defect_link_interior.append([i, j])
        else:
            if not is_edge_proper_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])
vts_defect_link = set_diff(unique(edges_defect_link), list_boundary)
count, count_interior = max(1, count), max(1, count_interior)

print("\nAll non-boundary edges: {}   -   reconstructed at a first try: {} %".
      format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
print("All interior edges: {}   -   reconstructed at a first try in the interior: {} %".
      format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))

log_file_detailed.write("\nNumber of non-boundary edges: {}, properly reconstructed at a first try: {} %".
                        format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file_detailed.write("\nNumber of interior edges: {}, properly reconstructed at a first try: {} %".
                        format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
log_file_detailed.flush()
log_file.write("\nNumber of non-boundary edges: {}, properly reconstructed at a first try: {} %".
               format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file.write("\nNumber of interior edges: {}, properly reconstructed at a first try: {} %".
               format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
log_file.flush()
log_file_detailed.write("\n\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

start_time = time()
print("\n\nApply a post process step to try to repair the triangulation...")
log_file_detailed.write("\n\nApply a post process step to try to repair the triangulation...")
log_file_detailed.flush()
log_file.write("\n\nApply a post process step to try to repair the triangulation...")
log_file.flush()

for vertex in range(num_pores):
    if vertex in ignore:
        continue
    if is_vertex_proper_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
        continue
    _, triangulation_reconstructed, _ = is_triangulatable(vertex, Matrix_out, list_boundary, Centroids, changed_tunnels,
                                                          predetermined=triangulation_reconstructed, post_process=True)

edges_defect_link, edges_defect_link_interior = [], []
count, count_interior = 0, 0
for i in range(num_pores):
    Nei = Matrix_out[i]
    for j in [item for item in Nei if item > i]:
        if sublist([i, j], union(list_boundary, list_boundary_vox)):
            continue
        count += 1
        if not intersect([i, j], union(list_boundary, list_boundary_vox)):
            count_interior += 1
            if not is_edge_proper_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])
                edges_defect_link_interior.append([i, j])
        else:
            if not is_edge_proper_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])
vts_defect_link = set_diff(unique(edges_defect_link), list_boundary)
count, count_interior = max(1, count), max(1, count_interior)

print("\nAll non-boundary edges: {}   -   reconstructed with a post process: {} %".
      format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
print("All interior edges: {}   -   reconstructed with a post process in the interior: {} %".
      format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))

log_file_detailed.write("\nNumber of non-boundary edges: {}, properly reconstructed with a post process: {} %".
                        format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file_detailed.write("\nNumber of interior edges: {}, properly reconstructed with a post process: {} %".
                        format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
log_file_detailed.flush()
log_file.write("\nNumber of non-boundary edges: {}, properly reconstructed with a post process: {} %".
               format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file.write("\nNumber of interior edges: {}, properly reconstructed with a post process: {} %".
               format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
log_file.flush()
log_file_detailed.write("\n\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


start_time = time()
print("\n\nChange links of non-proper edges and detect football pores and multi-connected pores...")
log_file_detailed.write("\n\nChange links of non-proper edges and detect football pores and multi-connected pores...")
log_file_detailed.flush()
log_file.write("\n\nChange links of non-proper edges and detect football pores and multi-connected pores...")
log_file.flush()

multi_edges = []
caused_by_digons = []

K2, K4 = [], []
for i in range(num_pores):
    Nei_i = Matrix_out[i]
    for j in [item for item in Nei_i if item > i]:
        if not sublist([i, j], union(list_boundary, list_boundary_vox)):
            K2.append([i, j])
        Nei_i_j = intersect(Nei_i, Matrix_out[j])
        for k in [item for item in Nei_i_j if item > j]:
            Nei_i_j_k = intersect(Nei_i_j, Matrix_out[k])
            for l in [item for item in Nei_i_j_k if item > k]:
                K4.append([i, j, k, l])

empty = []
processed = []
count = 0
loop_over = K2

while loop_over:
    edge = loop_over[0]
    count += 1
    multi = False
    print("\rProcessing: {} out of {}".format(count, len(edges_defect_link)), end="")
    if is_edge_proper_in_triangulation(edge, triangulation_reconstructed):
        loop_over = loop_over[1:]
        continue
    lst = mutual_neighbors(edge, Matrix_out)
    if len(lst) < 3 or (len(lst) < 6 and not cycles(edges_of(lst, Matrix_out), minimum_allowed=3)):
        caused_by_digons.append(edge)
        loop_over = loop_over[1:]
        continue
    vts = mutual_neighbors(edge, Matrix_out)
    edges = edges_of(vts, Matrix_out)
    Matrix_temp = remove_edge(edge[0], edge[1], Matrix_out)
    Holes_1 = holes(edge[0], Matrix_temp, list_boundary, Centroids, True)[0]
    Holes_1 = [hole for hole in Holes_1 if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
    Matrix_reduced = reduced_vertexlink(edge[0], Matrix_temp, list_boundary, Centroids, return_Matrix=True)[1]
    Holes_1_red = holes(edge[0], Matrix_reduced, list_boundary, Centroids, True)[0]
    Holes_1_red = [hole for hole in Holes_1_red if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
    Holes_1 += [hole for hole in Holes_1_red if not is_line_exact(hole, Holes_1)[0]]

    Holes_2 = holes(edge[1], Matrix_temp, list_boundary, Centroids, True)[0]
    Holes_2 = [hole for hole in Holes_2 if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
    Matrix_reduced = reduced_vertexlink(edge[1], Matrix_temp, list_boundary, Centroids, return_Matrix=True)[1]
    Holes_2_red = holes(edge[1], Matrix_reduced, list_boundary, Centroids, True)[0]
    Holes_2_red = [hole for hole in Holes_2_red if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
    Holes_2 += [hole for hole in Holes_2_red if not is_line_exact(hole, Holes_2)[0]]
    if len(Holes_1) > 1 and len(Holes_2) > 1:
        multi = True  # a probable multi edge
        if len(edges) > 25:  # too large to process
            if not is_line(edge, multi_edges):
                multi_edges.append(edge)
            loop_over = loop_over[1:]
            continue
        # try to fix the triangulation
        edges = edges_of(union(*(Holes_1 + Holes_2)), Matrix_out)
        C = cycles(edges, minimum_allowed=4, order_matters=True)
        edges = edges_of(union(*(Holes_1 + Holes_2), edge), Matrix_out)
        lst_remove = []
        for i in range(len(C)):
            if i in lst_remove:
                continue
            for j in range(i + 1, len(C)):
                if len(C[i]) != len(C[j]):
                    continue
                if set(C[j]) == set(C[i]) and is_line_EXACT(C[i], C[j]):
                    lst_remove.append(j)
        C = remove_index(C, lst_remove)
        lengths = unique([len(item) for item in C])
        C_ordered = []
        for i in lengths:
            C_ordered.extend([item for item in C if len(item) == i])
        candidates = []
        for cycle in C_ordered:
            Com_1 = components(edge[0], cycle, Matrix_temp)
            Com_2 = components(edge[1], cycle, Matrix_temp)
            if len(Com_1) != 1 or len(Com_2) != 1:
                continue
            Matrix_temp_sub = loads(dumps(Matrix_out))
            for i in Matrix_temp_sub[edge[0]]:
                if i in cycle:
                    continue
                Matrix_temp_sub = remove_edge(i, edge[1], Matrix_temp_sub)
            Holes = holes(edge[0], Matrix_temp_sub, list_boundary, Centroids, True)[0]
            Holes = [hole for hole in Holes if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
            if Holes:
                continue
            Matrix_temp_sub = loads(dumps(Matrix_out))
            for i in Matrix_temp_sub[edge[1]]:
                if i in cycle:
                    continue
                Matrix_temp_sub = remove_edge(i, edge[0], Matrix_temp_sub)
            Holes = holes(edge[1], Matrix_temp_sub, list_boundary, Centroids, True)[0]
            Holes = [hole for hole in Holes if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
            if Holes:
                continue
            candidates.append(cycle)
        lst_remove = []
        for index in range(len(candidates)):
            if len(set_diff(range(len(candidates)), lst_remove)) <= 10:
                break
            cycle = candidates[index]
            d = sum([is_line(item, cycle) for item in K4])
            if sum([is_line(edge, cycle) for edge in edges]) - d >= 2 * len(cycle) - 3:  # hole is already triangulated
                lst_remove.append(index)
                continue
            Matrix_temp_sub = loads(dumps(Matrix_out))
            for i in Matrix_temp_sub[edge[0]]:
                if i in cycle:
                    continue
                Matrix_temp_sub = remove_edge(i, edge[1], Matrix_temp_sub)
            Centroids_temp = loads(dumps(Centroids))
            center = [0, 0, 0]
            for i in cycle:
                center += np.array(Centroids[i])
            center /= len(cycle)
            center = 2 * center - np.array(Centroids[edge[0]])
            Centroids_temp[edge[1]] = list(center)
            Holes = holes(edge[0], Matrix_temp_sub, list_boundary, Centroids_temp, True)[0]
            Holes = [hole for hole in Holes if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
            if Holes:
                lst_remove.append(index)
                continue
            Matrix_temp_sub = loads(dumps(Matrix_out))
            for i in Matrix_temp_sub[edge[1]]:
                if i in cycle:
                    continue
                Matrix_temp_sub = remove_edge(i, edge[0], Matrix_temp_sub)
            Centroids_temp = loads(dumps(Centroids))
            center = [0, 0, 0]
            for i in cycle:
                center += np.array(Centroids[i])
            center /= len(cycle)
            center = 2 * center - np.array(Centroids[edge[1]])
            Centroids_temp[edge[0]] = list(center)
            Holes = holes(edge[1], Matrix_temp_sub, list_boundary, Centroids_temp, True)[0]
            Holes = [hole for hole in Holes if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
            if Holes:
                lst_remove.append(index)
                continue
        candidates = remove_index(candidates, lst_remove)
        print("\nnumber of candidates for {} = {}".format(edge, len(candidates)))
        if not candidates or len(candidates) > 20:
            if multi and not is_line(edge, multi_edges):
                multi_edges.append(edge)
            loop_over = loop_over[1:]
            continue
    else:
        if len(edges) > 25:  # too large to process
            loop_over = loop_over[1:]
            continue
        if Holes_1 + Holes_2:
            C = cycles(edges, minimum_allowed=4, order_matters=True)
        else:
            C = cycles(edges, minimum_allowed=3, order_matters=True)
        edges = edges_of(union(*(Holes_1 + Holes_2), edge), Matrix_out)
        lengths = unique([len(item) for item in C])
        C_ordered = []
        for i in lengths:
            C_ordered.extend([item for item in C if len(item) == i])
        candidates = []
        for cycle in C_ordered:
            Com_1 = components(edge[0], cycle, Matrix_temp)
            Com_2 = components(edge[1], cycle, Matrix_temp)
            if len(Com_1) != 1 or len(Com_2) != 1:
                continue
            Matrix_temp_sub = loads(dumps(Matrix_out))
            for i in Matrix_temp_sub[edge[0]]:
                if i in cycle:
                    continue
                Matrix_temp_sub = remove_edge(i, edge[1], Matrix_temp_sub)
            Holes = holes(edge[0], Matrix_temp_sub, list_boundary, Centroids, True)[0]
            Holes = [hole for hole in Holes if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
            if Holes:
                continue
            Matrix_temp_sub = loads(dumps(Matrix_out))
            for i in Matrix_temp_sub[edge[1]]:
                if i in cycle:
                    continue
                Matrix_temp_sub = remove_edge(i, edge[0], Matrix_temp_sub)
            Holes = holes(edge[1], Matrix_temp_sub, list_boundary, Centroids, True)[0]
            Holes = [hole for hole in Holes if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
            if Holes:
                continue
            candidates.append(cycle)
        lst_remove = []
        for index in range(len(candidates)):
            if len(set_diff(range(len(candidates)), lst_remove)) <= 10:
                break
            cycle = candidates[index]
            d = sum([is_line(item, cycle) for item in K4])
            if len(cycle) > 3 and sum([is_line(edge, cycle) for edge in edges]) - d >= 2 * len(cycle) - 3:  # hole is already triangulated
                lst_remove.append(index)
            Matrix_temp_sub = loads(dumps(Matrix_out))
            for i in Matrix_temp_sub[edge[0]]:
                if i in cycle:
                    continue
                Matrix_temp_sub = remove_edge(i, edge[1], Matrix_temp_sub)
            Centroids_temp = loads(dumps(Centroids))
            center = [0, 0, 0]
            for i in cycle:
                center += np.array(Centroids[i])
            center /= len(cycle)
            center = 2 * center - np.array(Centroids[edge[0]])
            Centroids_temp[edge[1]] = list(center)
            Holes = holes(edge[0], Matrix_temp_sub, list_boundary, Centroids_temp, True)[0]
            Holes = [hole for hole in Holes if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
            if Holes:
                lst_remove.append(index)
                continue
            Matrix_temp_sub = loads(dumps(Matrix_out))
            for i in Matrix_temp_sub[edge[1]]:
                if i in cycle:
                    continue
                Matrix_temp_sub = remove_edge(i, edge[0], Matrix_temp_sub)
            Centroids_temp = loads(dumps(Centroids))
            center = [0, 0, 0]
            for i in cycle:
                center += np.array(Centroids[i])
            center /= len(cycle)
            center = 2 * center - np.array(Centroids[edge[1]])
            Centroids_temp[edge[0]] = list(center)
            Holes = holes(edge[1], Matrix_temp_sub, list_boundary, Centroids_temp, True)[0]
            Holes = [hole for hole in Holes if all([is_line(e, edges) for e in edges_of_cycle(hole)])]
            if Holes:
                lst_remove.append(index)
                continue
        candidates = remove_index(candidates, lst_remove)
        print("\nnumber of candidates for {} = {}".format(edge, len(candidates)))
        if not candidates or len(candidates) > 20:
            candidates = Holes_1 + Holes_2
            if candidates:
                print("number of candidates reduced to", len(Holes_1 + Holes_2))
        if not candidates:
            loop_over = loop_over[1:]
            continue
    # we consider the second neighbors
    lst_vertices = union(Matrix_out[edge[0]], Matrix_out[edge[1]])
    lst_vertices = unique([Matrix_out[item] for item in lst_vertices])
    lst_edges = edges_of(lst_vertices, Matrix_out)
    proper_edges_initial = [edge for edge in lst_edges if is_edge_proper_in_triangulation(edge,
                                                                                          triangulation_reconstructed)]
    proper_vertices_initial = [v for v in lst_vertices if is_vertex_proper_in_triangulation(
        v, triangulation_reconstructed, Matrix_out[v])]
    proper_vertices_initial_int = set_diff(proper_vertices_initial, list_boundary)
    gains_edges, gains_vertices, gains_vertices_int, retriangulations, lsts, candidates_unfold, emptys, losses = [], [], [], [], [], [], [], []
    for cycle in candidates:
        triangulation_test = loads(dumps(triangulation_reconstructed))
        Ind = []
        for u in mutual_neighbors(edge, Matrix_out):
            Ind = union(Ind, ind([edge[0], u], triangulation_test))
            Ind = union(Ind, ind([edge[1], u], triangulation_test))
        for edge_sub in combinations(cycle, 2):
            Ind = union(Ind, ind(edge_sub, triangulation_test))
        triangulation_test = remove_index(triangulation_test, Ind)
        Centroids_temp_0, Centroids_temp_1 = loads(dumps(Centroids)), loads(dumps(Centroids))
        center = [0, 0, 0]
        for i in cycle:
            center += np.array(Centroids[i])
        center /= len(cycle)
        center_0 = 2 * center - np.array(Centroids[edge[1]])
        center_1 = 2 * center - np.array(Centroids[edge[0]])
        Centroids_temp_0[edge[1]] = list(center_1)
        Centroids_temp_1[edge[0]] = list(center_0)

        edges = edges_of_cycle(cycle)
        star = [edge + item for item in edges]
        empty_sub = [K4[item] for item in ind(edge, K4) if not is_line(K4[item], star)]
        empty_test = empty + [item for item in empty_sub if not is_line(item, empty)]
        triangles = []
        if len(cycle) > 3:
            triangles = cycles_deg(edges_of(cycle, Matrix_out), 3)
        if triangles and len(triangles) < 4:
            for num in range(2 ** len(triangles)):
                candidates_unfold.append(cycle)
                pos = bin(num)[2:].zfill(len(triangles))
                empty_test_sub = loads(dumps(empty_test))
                retriangulation = triangulation_test + star
                for index_sub in range(len(pos)):
                    if pos[index_sub] == "1":
                        empty_test_sub += [item for item in [[edge[0]] + triangles[index_sub],
                                           [edge[1]] + triangles[index_sub]] if not is_line(item, empty_test_sub)]
                emptys.append(empty_test_sub)
                for v_sub in union(unique(candidates), edge):
                    retriangulation = is_triangulatable(v_sub, Matrix_out, list_boundary, Centroids,
                                                        changed_tunnels, predetermined=retriangulation,
                                                        predetermined_empty=empty_test_sub)[1]
                retriangulations.append(retriangulation)

                proper_vertices_final = [v for v in lst_vertices if is_vertex_proper_in_triangulation(
                    v, retriangulation, Matrix_out[v])]
                proper_edges_final = [edge for edge in lst_edges if is_edge_proper_in_triangulation(edge,
                                                                                                    retriangulation)]
                proper_vertices_final_int = set_diff(proper_vertices_final, list_boundary)
                lsts.append(proper_edges_final)
                losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
                gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
                gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
                gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))
        else:
            candidates_unfold.append(cycle)
            emptys.append(empty_test)
            retriangulation = triangulation_test + star
            for v_sub in union(unique(candidates), edge):
                retriangulation = is_triangulatable(v_sub, Matrix_out, list_boundary, Centroids,
                                                    changed_tunnels, predetermined=retriangulation,
                                                    predetermined_empty=empty_test)[1]
            retriangulations.append(retriangulation)

            proper_vertices_final = [v for v in lst_vertices if is_vertex_proper_in_triangulation(
                v, retriangulation, Matrix_out[v])]
            proper_edges_final = [edge for edge in lst_edges if is_edge_proper_in_triangulation(edge, retriangulation)]
            proper_vertices_final_int = set_diff(proper_vertices_final, list_boundary)
            lsts.append(proper_edges_final)
            losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
            gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
            gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
            gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))
        candidates_unfold.append(cycle)
        emptys.append(empty_test)
        retriangulation = triangulation_test + star
        for v_sub in union(unique(candidates), edge):
            Centroid_current = Centroids
            if v_sub == edge[0]:
                Centroid_current = Centroids_temp_0
            elif v_sub == edge[1]:
                Centroid_current = Centroids_temp_1
            retriangulation = is_triangulatable(v_sub, Matrix_out, list_boundary, Centroid_current,
                                                changed_tunnels, predetermined=retriangulation,
                                                predetermined_empty=empty_test)[1]
        retriangulations.append(retriangulation)

        proper_vertices_final = [v for v in lst_vertices if is_vertex_proper_in_triangulation(v, retriangulation,
                                                                                              Matrix_out[v])]
        proper_edges_final = [edge for edge in lst_edges if is_edge_proper_in_triangulation(edge, retriangulation)]
        proper_vertices_final_int = set_diff(proper_vertices_final, list_boundary)
        lsts.append(proper_edges_final)
        losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
        gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
        gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
        gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))

    m = max(gains_edges)
    if m > 0:
        Ind = [index for index in range(len(gains_edges)) if gains_edges[index] == m]
        print("gains_edges: {}  -  {}".format(m, Ind))
        for index in set_diff(range(len(gains_vertices)), Ind):
            gains_vertices[index] = -math.inf
        m = max(gains_vertices)
        Ind = [index for index in range(len(gains_vertices)) if gains_vertices[index] == m]
        print("gains_vertices: {}  -  {}".format(m, Ind))
        for index in set_diff(range(len(gains_vertices_int)), Ind):
            gains_vertices_int[index] = -math.inf
        Ind = gains_vertices_int.index(max(gains_vertices_int))
        print("gains_vertices_interior: {}  -  {}".format(m, [Ind]))
        triangulation_reconstructed = loads(dumps(retriangulations[Ind]))
        empty = loads(dumps(emptys[Ind]))
        loop_over.extend(losses[Ind])
        print("best found:", candidates_unfold[Ind])
        log_file_detailed.write("\nConsidered edge: {}, best link found: {}, number of edges gained: {}"
                                .format(edge, candidates_unfold[Ind], max(gains_edges)))
        log_file_detailed.flush()
    else:
        if multi and not is_line(edge, multi_edges):
            multi_edges.append(edge)
    loop_over = loop_over[1:]


edges_defect_link, edges_defect_link_interior = [], []
count, count_interior, count_reconstructable, count_reconstructable_interior = 0, 0, 0, 0
for i in range(num_pores):
    Nei = Matrix_out[i]
    for j in [item for item in Nei if item > i]:
        if sublist([i, j], union(list_boundary, list_boundary_vox)):
            continue
        count += 1
        if not intersect([i, j], union(list_boundary, list_boundary_vox)):
            count_interior += 1
            if not is_line([i, j], multi_edges + caused_by_digons):
                count_reconstructable += 1
                count_reconstructable_interior += 1
            if not is_edge_proper_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])
                edges_defect_link_interior.append([i, j])
        else:
            if not is_line([i, j], multi_edges + caused_by_digons):
                count_reconstructable += 1
            if not is_edge_proper_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])
count, count_interior = max(1, count), max(1, count_interior)
count_reconstructable, count_reconstructable_interior = max(1, count_reconstructable), max(1, count_reconstructable_interior)


print("\nAll non-boundary edges: {}   -   reconstructed after changing links of non-proper edges: {} %".
      format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
print("All interior edges: {}   -   reconstructed after changing links of non-proper edges (interior): {} %".
      format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
print("Number of multi-edges:", len(multi_edges))
print("Number of non-properly reconstructed edges around football pores:", len(caused_by_digons))
leftovers = [item for item in edges_defect_link if not is_line(item, multi_edges + caused_by_digons)]
leftovers_interior = [item for item in edges_defect_link_interior if not is_line(item, multi_edges + caused_by_digons)]
print("All reconstructable non-boundary edges: {}   -   reconstructed among reconstructable edges: {} %".
      format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
print("All reconstructable interior edges: {}   -   reconstructed among reconstructable edges (interior): {} %".
      format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))


log_file_detailed.write("\n\nNumber of non-boundary edges: {}, properly reconstructed after changing links of non-proper edges: {} %".
                        format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file_detailed.write("\nNumber of interior edges: {}, properly reconstructed after changing links of non-proper edges: {} %".
                        format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
log_file.write("\nNumber of non-boundary edges: {}, properly reconstructed after changing links of non-proper edges: {} %".
               format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file.write("\nNumber of interior edges: {}, properly reconstructed after changing links of non-proper edges: {} %".
               format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
s = "pairs"
if len(multi_edges) == 1:
    s = "pair"
log_file_detailed.write("\n\nNumber of multi-connected pores: {} {}".format(len(multi_edges), s))
log_file_detailed.write("\nNumber of non-properly reconstructed edges around football pores: {}".format(len(caused_by_digons)))
log_file_detailed.write("\n\nNumber of reconstructable non-boundary edges: {}, properly reconstructed among reconstructable edges: {} %".
                        format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
log_file_detailed.write("\nNumber of reconstructable interior edges: {}, properly reconstructed after changing links of non-proper edges: {} %".
                        format(count_reconstructable_interior,
                               round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))
log_file_detailed.write("\n\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nNumber of multi-connected pores: {} {}".format(len(multi_edges), s))
log_file.write("\nNumber of non-properly reconstructed edges around football pores: {}".format(len(caused_by_digons)))
log_file.write("\nNumber of reconstructable non-boundary edges: {}, properly reconstructed among reconstructable edges: {} %".
               format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
log_file.write("\nNumber of reconstructable interior edges: {}, properly reconstructed after changing links of non-proper edges: {} %".
               format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


start_time = time()
print("\n\nChange tunnels in the triangulation...")
log_file_detailed.write("\n\nChange tunnels in the triangulation...")
log_file_detailed.flush()
log_file.write("\n\nChange tunnels in the triangulation...")
log_file.flush()

K5 = find_K5s(Matrix_out)

loop_over = [v for v in ordered_pores if is_line(v, leftovers)]
while loop_over:
    vertex = loop_over[0]
    T, problems = is_vertex_proper_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex],
                                                    return_problems=True)
    if T:
        loop_over = loop_over[1:]
        continue
    problems = [item for item in problems if not is_line([vertex, item], multi_edges + caused_by_digons)]
    if not problems:
        loop_over = loop_over[1:]
        continue
    updated = False
    K5_all = [K5[item] for item in ind(vertex, K5)]
    lst_vertices = unique(K5_all)
    proper_vertices_initial = [v for v in lst_vertices if is_vertex_proper_in_triangulation(v, triangulation_reconstructed,
                                                                                            Matrix_out[v])]
    proper_vertices_initial_int = set_diff(proper_vertices_initial, list_boundary)
    K5_targeted = [item for item in K5_all if intersect(item, problems)]
    for K5_current in K5_targeted:
        if is_vertex_proper_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
            continue
        lst_edges, tets = [], []
        for clique in K5_all:
            if len(intersect(clique, K5_current)) >= 4:
                lst_edges += [list(item) for item in combinations(clique, 2) if not is_line(item, lst_edges)]
                tets += [list(item) for item in combinations(clique, 4) if not is_line(item, tets)]
        proper_edges_initial = [edge for edge in lst_edges if is_edge_proper_in_triangulation(edge,
                                                                                           triangulation_reconstructed)]
        gains_edges, gains_vertices, gains_vertices_int, retriangulations, lsts, emptys, losses = [], [], [], [], [], [], []
        triangulation_test = loads(dumps(triangulation_reconstructed))
        tets = [list(item) for item in combinations(K5_current, 4)]
        for tet in tets:
            triangulation_test = remove_index(triangulation_test, ind(tet, triangulation_test))
        for num in range(2 ** 5 - 1):  # K5 with 'all-ones' is not possible, so we disregard 31 = 11111
            pos = bin(num)[2:].zfill(5)
            retriangulation = loads(dumps(triangulation_test))
            retriangulation += [tets[index] for index in range(len(pos)) if pos[index] == "1"]
            emp = [tets[index] for index in range(len(pos)) if pos[index] == "0"]
            empty_test = empty + [item for item in emp if not is_line(item, empty)]
            emptys.append(loads(dumps(empty_test)))
            retriangulations.append(retriangulation)

            proper_vertices_final = [v for v in lst_vertices if is_vertex_proper_in_triangulation(v, retriangulation,
                                                                                                  Matrix_out[v])]
            proper_edges_final = [edge for edge in lst_edges if is_edge_proper_in_triangulation(edge, retriangulation)]
            proper_vertices_final_int = set_diff(proper_vertices_final, list_boundary)
            lsts.append(proper_edges_final)
            losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
            gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
            gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
            gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))

        m = max(gains_edges)
        if m > 0:
            print("vertex: {}, 5-clique: {}".format(vertex, K5_current))
            Ind = [index for index in range(len(gains_edges)) if gains_edges[index] == m]
            print("gains_edges: {}  -  {}".format(m, Ind))
            for index in set_diff(range(len(gains_vertices)), Ind):
                gains_vertices[index] = -math.inf
            m = max(gains_vertices)
            Ind = [index for index in range(len(gains_vertices)) if gains_vertices[index] == m]
            print("gains_vertices: {}  -  {}".format(m, Ind))
            for index in set_diff(range(len(gains_vertices_int)), Ind):
                gains_vertices_int[index] = -math.inf
            Ind = gains_vertices_int.index(max(gains_vertices_int))
            print("gains_vertices_interior: {}  -  {}".format(m, [Ind]))
            triangulation_reconstructed = loads(dumps(retriangulations[Ind]))
            empty = loads(dumps(emptys[Ind]))
            loop_over.extend(unique(losses[Ind]))
            updated = True
            log_file_detailed.write("\nConsidered 5-clique: {}, number of edges gained after a tunnel change: {}".
                                    format(K5_current, max(gains_edges)))
            log_file_detailed.flush()
    if not updated:
        loop_over = loop_over[1:]

edges_defect_link, edges_defect_link_interior = [], []
count, count_interior, count_reconstructable, count_reconstructable_interior = 0, 0, 0, 0
for i in range(num_pores):
    Nei = Matrix_out[i]
    for j in [item for item in Nei if item > i]:
        if sublist([i, j], union(list_boundary, list_boundary_vox)):
            continue
        count += 1
        if not intersect([i, j], union(list_boundary, list_boundary_vox)):
            count_interior += 1
            if not is_line([i, j], multi_edges + caused_by_digons):
                count_reconstructable += 1
                count_reconstructable_interior += 1
            if not is_edge_proper_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])
                edges_defect_link_interior.append([i, j])
        else:
            if not is_line([i, j], multi_edges + caused_by_digons):
                count_reconstructable += 1
            if not is_edge_proper_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])

count, count_interior = max(1, count), max(1, count_interior)
count_reconstructable, count_reconstructable_interior = max(1, count_reconstructable), max(1, count_reconstructable_interior)

print("\nAll non-boundary edges: {}   -   reconstructed after changing tunnels: {} %".
      format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
print("All interior edges: {}   -   reconstructed after changing tunnels (interior): {} %".
      format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
print("Number of multi-edges:", len(multi_edges))
print("Number of non-properly reconstructed edges around football pores:", len(caused_by_digons))
leftovers = [item for item in edges_defect_link if not is_line(item, multi_edges + caused_by_digons)]
leftovers_interior = [item for item in edges_defect_link_interior if not is_line(item, multi_edges + caused_by_digons)]
print("All reconstructable non-boundary edges: {}   -   reconstructed among reconstructable edges: {} %".
      format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
print("All reconstructable interior edges: {}   -   reconstructed among reconstructable edges (interior): {} %".
      format(count_reconstructable_interior,
             round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))


log_file_detailed.write("\n\nNumber of non-boundary edges: {}, properly reconstructed after changing tunnels: {} %".
                        format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file_detailed.write("\nNumber of interior edges: {}, properly reconstructed after changing tunnels: {} %".
                        format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
log_file.write("\nNumber of non-boundary edges: {}, properly reconstructed after changing tunnels: {} %".
               format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file.write("\nNumber of interior edges: {}, properly reconstructed after changing tunnels: {} %".
               format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))

log_file_detailed.write("\n\nNumber of reconstructable non-boundary edges: {}, properly reconstructed among reconstructable edges: {} %".
                        format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
log_file_detailed.write("\nNumber of reconstructable interior edges: {}, properly reconstructed among reconstructable edges: {} %".
                        format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))

log_file.write("\nNumber of reconstructable non-boundary edges: {}, properly reconstructed among reconstructable edges: {} %".
               format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
log_file.write("\nNumber of reconstructable interior edges: {}, properly reconstructed among reconstructable edges: {} %".
               format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))


start_time = time()
print("\n\nTry to remove excessive edges to repair defect links...")
log_file_detailed.write("\n\nTry to remove excessive edges to repair defect links...")
log_file_detailed.flush()
log_file.write("\n\nTry to remove excessive edges to repair defect links...")
log_file.flush()

vertices_defect_link = []
for vertex in range(num_pores):
    if vertex in union(list_boundary, list_boundary_vox, ignore):
        continue
    if vertex in union(unique(caused_by_digons), unique(multi_edges)):
        continue
    if not is_vertex_proper_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
        vertices_defect_link.append(vertex)

for vertex in vertices_defect_link:
    if is_vertex_proper_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
        continue
    candidates, improvements, distances = [], [], []
    for u in Matrix_out[vertex]:
        Matrix_temp = remove_edge(vertex, u, Matrix_out)
        triangulation_temp = remove_index(triangulation_reconstructed, ind([vertex, u], triangulation_reconstructed))
        if not is_vertex_proper_in_triangulation(vertex, triangulation_temp, Matrix_temp[vertex]):
            continue
        lst = union(mutual_neighbors([vertex, u], Matrix_out), u)
        T1 = [is_vertex_proper_in_triangulation(v, triangulation_reconstructed, Matrix_out[v]) for v in lst]
        T2 = [is_vertex_proper_in_triangulation(v, triangulation_temp, Matrix_temp[v]) for v in lst]
        if not all([not T1[item] or T2[item] for item in range(len(T1))]):
            continue
        candidates.append(u)
        improvements.append(sum(T2) - sum(T1))
        distances.append(distance_vox(vertex, u, Centroids, coordinates, voxels, dist_path))
    if not candidates:
        continue
    Ind = [item for item in range(len(candidates)) if improvements[item] == max(improvements)]
    Ind = distances.index(max([distances[item] for item in Ind]))
    u, dist = candidates[Ind], distances[Ind]
    Matrix_out = remove_edge(vertex, u, Matrix_out)
    Register_Removed.append([vertex, u])
    log_file_detailed.write(
        "\nRemoved excessive edge {} of length {} to repair the link of {}".format([vertex, u], dist, vertex))
    log_file_detailed.flush()
    print("Removed excessive edge {} of length {} to repair the link of {}".format([vertex, u], dist, vertex))
    triangulation_reconstructed = remove_index(triangulation_reconstructed,
                                               ind([vertex, u], triangulation_reconstructed))

vertices_defect_link = []
lst_good = []
count_vts = 0
for vertex in range(num_pores):
    if vertex in union(list_boundary, list_boundary_vox, ignore):
        continue
    if vertex in union(unique(caused_by_digons), unique(multi_edges)):
        continue
    count_vts += 1
    if is_vertex_proper_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
        lst_good.append(vertex)
    else:
        vertices_defect_link.append(vertex)

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
avg_deg += sum([len(item) for item in football_pores_neighbors])

edges_defect_link, edges_defect_link_interior = [], []
count, count_interior = 0, 0
for i in range(num_pores):
    Nei = Matrix_out[i]
    for j in [item for item in Nei if item > i]:
        if sublist([i, j], union(list_boundary, list_boundary_vox)):
            continue
        if intersect([i, j], union(unique(caused_by_digons), unique(multi_edges))):
            continue
        count += 1
        if not intersect([i, j], union(list_boundary, list_boundary_vox)):
            count_interior += 1
            if not is_edge_proper_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])
                edges_defect_link_interior.append([i, j])
        else:
            if not is_edge_proper_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])
vts_defect_link = set_diff(unique(edges_defect_link), list_boundary)
count, count_interior = max(1, count), max(1, count_interior)
##
print("Num tetrahedra:", len(triangulation_reconstructed))
##
print("Multi-connected pairs of pores: {}".format(len(multi_edges)))
print("Number of non-properly reconstructed edges around football pores: {}".format(len(caused_by_digons)))
print("Non-boundary edges: {}".format(count))
print("Edges in the strict interior: {}".format(count_interior))
print("Percentage of properly reconstructed edges (non-boundary) {} %".format(
    round(100 * (1 - len(edges_defect_link) / count), 2)))
print("Percentage of properly reconstructed edges (strict interior) {} %".format(
    round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
print("\nAll pores:", num_pores)
print("All interior vertices:", len(set_diff(range(num_pores), union(list_boundary, list_boundary_vox))))
print("Reconstructable interior vertices (among reconstructables): {} %".format(
    round(100 * (1 - len(vertices_defect_link) / count_vts), 2)))
print("Average face degree of interior pores = {}".format(round(avg_deg / len(list_interior), 3)))

log_file_detailed.write("\n\nMulti-connected pairs of pores: {}".format(len(multi_edges)))
log_file_detailed.write("\nNumber of non-properly reconstructed edges around football pores: {}".
                        format(len(caused_by_digons)))
log_file_detailed.write("\nNon-boundary edges: {}".format(count))
log_file_detailed.write("\nEdges in the strict interior: {}".format(count_interior))
log_file_detailed.write("\nPercentage of properly reconstructed edges (non-boundary) {} %".format(
    round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file_detailed.write("\nPercentage of properly reconstructed edges (strict interior) {} %".format(
    round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
log_file_detailed.write("\n\nAll pores: {}".format(num_pores))
log_file_detailed.write("\nAll interior vertices: {}".format(len(set_diff(range(num_pores), union(list_boundary,
                                                                                                  list_boundary_vox)))))
log_file_detailed.write("\nReconstructable interior vertices (among reconstructables): {} %".format(
    round(100 * (1 - len(vertices_defect_link) / count_vts), 2)))
log_file_detailed.write("\nAverage face degree of interior pores = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n\n{}\n".format(compute_time(start_time, True)))
log_file_detailed.flush()

log_file.write("\n\nMulti-connected pairs of pores: {}".format(len(multi_edges)))
log_file.write("\nNumber of non-properly reconstructed edges around football pores: {}".format(len(caused_by_digons)))
log_file.write("\nNon-boundary edges: {}".format(count))
log_file.write("\nEdges in the strict interior: {}".format(count_interior))
log_file.write("\nPercentage of properly reconstructed edges (non-boundary) {} %".format(
    round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file.write("\nPercentage of properly reconstructed edges (strict interior) {} %".format(
    round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
log_file.write("\n\nAll pores: {}".format(num_pores))
log_file.write("\nAll interior vertices: {}".format(len(set_diff(range(num_pores),
                                                                 union(list_boundary, list_boundary_vox)))))
log_file.write("\nReconstructable interior vertices (among reconstructables): {} %".format(
    round(100 * (1 - len(vertices_defect_link) / count_vts), 2)))
log_file.write("\nAverage face degree of interior pores = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n\n{}\n".format(compute_time(start_time, True)))
log_file.flush()

degrees = {}
lst = union(list_boundary, list_boundary_vox)
for vertex in range(len(Matrix_out)):
    if vertex in lst:
        continue
    if len(Matrix_out[vertex]) in degrees:
        degrees[len(Matrix_out[vertex])] += 1
    else:
        degrees[len(Matrix_out[vertex])] = 1
count = 0
for i in range(4, 400):
    if i in degrees:
        count += degrees[i]
for i in range(4, 400):
    if i in degrees:
        print("{} vertices: {}%".format(i, round(100 * degrees[i] / count, 2)))
        log_file_detailed.write("\n{} vertices: {}%".format(i, round(100 * degrees[i] / count, 2)))
        log_file_detailed.flush()
        log_file.write("\n{} vertices: {}%".format(i, round(100 * degrees[i] / count, 2)))
        log_file.flush()

save_list_of_lists(Matrix_out, "adjacency_matrix.txt")
save_list(football_pores, "football_pores.txt")
save_list_of_lists(football_pores_neighbors, "football_pores_neighbors.txt")
save_list_of_lists(changed_tunnels, "changed_tunnels.txt")
save_list_of_lists(caused_by_digons, "caused_by_football_pores.txt")
save_list_of_lists(multi_edges, "multi_edges.txt")
save_list_of_lists(triangulation_reconstructed, "triangulation_reconstructed.txt")

log_file_detailed.close()
log_file.close()
