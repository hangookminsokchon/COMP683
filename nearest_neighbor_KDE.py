'''
script to compare distances between two groups of points. Specifically, for each
point in group 1, finds the point in group 2 that is the smallest distance away
and outputs the distances and names of the points in group 2 that meet these
criteria for the points in group 1

This version assumes a csv file called 'all_cells.csv'
'''

import numpy as np
import pandas as pd
   
# ADDITIONAL
#--------------------------------------------------------------------------------------------------------------
def build_cell_class_dict(coord_dict, class_map):
    grid_res = class_map.shape[0]  # assume square

    xs = [coord[0] for coord in coord_dict.values()]
    ys = [coord[1] for coord in coord_dict.values()]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    cell_class_dict = {}
    for cell, (xv, yv) in coord_dict.items():
        # normalize to [0,1]
        tx = (xv - xmin) / (xmax - xmin) if xmax > xmin else 0.0
        ty = (yv - ymin) / (ymax - ymin) if ymax > ymin else 0.0
        # clamp
        tx = min(max(tx, 0.0), 1.0)
        ty = min(max(ty, 0.0), 1.0)

        # grid index (nearest)
        j = int(round(tx * (grid_res - 1)))  # col
        i = int(round(ty * (grid_res - 1)))  # row

        cell_class_dict[cell] = int(class_map[i, j])

    return cell_class_dict


def compute_distances_KDE(cluster_dict, coord_dict, cluster1, cluster2, cell_class_dict):
    distance_dict = {}
    print(len(cluster_dict[cluster1]))
    print(len(cluster_dict[cluster2]))

    for cell_number, cell1 in enumerate(cluster_dict[cluster1]):
        distance_dict[cell1] = ['nonsense', float('Inf')]

        if cell_number % 100 == 0:
            print(cell_number / len(cluster_dict[cluster1]))

        class1 = cell_class_dict[cell1]

        for cell2 in cluster_dict[cluster2]:
            if cell_class_dict[cell2] == class1:
                continue

            x_distance = abs(coord_dict[cell1][0] - coord_dict[cell2][0])
            y_distance = abs(coord_dict[cell1][1] - coord_dict[cell2][1])
            total_distance = (x_distance**2 + y_distance**2)**0.5

            if total_distance < distance_dict[cell1][1]:
                distance_dict[cell1] = [cell2, total_distance]

    return distance_dict

# Load computed classmap
class_map = pd.read_csv('all_cells_example_classmap.csv')
class_map = np.array(class_map)

coord_dict, cluster_dict = get_coords('all_cells_example2.csv')
valid_clusters = sorted(list(cluster_dict.keys()))
cell_class_dict = build_cell_class_dict(coord_dict, class_map)

for first_cluster in valid_clusters:
    for second_cluster in valid_clusters:
        if first_cluster != second_cluster:
            print(first_cluster, second_cluster)
            distance_dict = compute_distances_KDE(
                cluster_dict, coord_dict, first_cluster, second_cluster, cell_class_dict
            )
            output_distances(distance_dict, first_cluster, second_cluster) 
   
