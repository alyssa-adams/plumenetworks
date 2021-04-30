import csv
import re
import scipy
import copy
import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite

from sklearn import datasets, linear_model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt

from pybdm import BDM
import pybdm
pybdm.options.set(raise_if_zero=False)


if __name__ == '__main__':

    # -------------- Load in data ----------------

    # read in network plume file
    reaction_data = {}

    with open('data/Plume_chemical_reaction_table.xlsx - Sheet1.csv') as csvfile:
        reader = csv.reader(csvfile)

        # skip headers
        next(reader)
        next(reader)

        for row in reader:

            metabolism = row[0]
            reaction_id = row[1]
            if len(reaction_id) == 0:
                continue
            reaction = row[3]
            electrons = int(row[4])
            dg_1 = float(row[5])
            dg_25 = float(row[6])
            dg_100 = float(row[7])

            ins = reaction.split('â†’')[0]
            outs = reaction.split('â†’')[1]
            ins = ins.split(' + ')
            outs = outs.split(' + ')
            ins = list(map(lambda x: x.strip(), ins))
            outs = list(map(lambda x: x.strip(), outs))

            reaction_data[reaction_id] = {}
            reaction_data[reaction_id]['metabolism'] = metabolism
            reaction_data[reaction_id]['electrons'] = electrons
            reaction_data[reaction_id]['dg_1'] = dg_1
            reaction_data[reaction_id]['dg_25'] = dg_25
            reaction_data[reaction_id]['dg_100'] = dg_100
            reaction_data[reaction_id]['ins'] = ins
            reaction_data[reaction_id]['outs'] = outs

    # load in the different plume chemical network data
    plumes = ['guaymas_basin', 'cayman_von_damm', 'cayman_picard', 'lau_basin_abe', 'lau_basin_kilo_moana', 'lau_basin_mariner']
    temps = ['4.9C', '2.95C']

    # enter the node attribute data for each
    plume_data = {}
    plume_data_percent_yield = {}

    plume_data['guaymas_basin'] = {
        '4.9C': {
            'total_energy': 20.79,
            'R-H-1': 4.399017142,
            'R-M-1': 0.4343748767,
            'R-S-1': 0.00,
            'R-S-2': 15.91578707,
            'R-I-2': 0.03643534153,
            'N/A': 0.00
                },
        '2.95C': {
            'total_energy': 3.28,
            'R-H-1': 0.04398480002,
            'R-M-1': 3.026942646,
            'R-S-1': 0.2098135263,
            'R-S-2': 0.00,
            'R-I-2': 0.00,
            'N/A': 0.00
                },
        }

    plume_data['cayman_von_damm'] = {
        '4.9C': {
            'total_energy': 104.65,
            'R-H-1': 52.95336756,
            'R-M-1': 27.06828318,
            'R-S-1': 17.00712931,
            'R-S-2': 6.173328209,
            'R-A-1': 0.632495608,
            'N/A': 0.818491034
        },
        '2.95C': {
            'total_energy': 33.58,
            'R-H-1': 17.16413769,
            'R-M-1': 8.867828407,
            'R-S-1': 7.331612582,
            'R-S-2': 0.00,
            'R-A-1': 0.218907133,
            'N/A': 0.00
        },
    }

    plume_data['cayman_picard'] = {
        '4.9C': {
            'total_energy': 82.31,
            'R-H-1': 32.1528303,
            'R-M-1': 0.657333352,
            'R-S-1': 47.39114927,
            'R-A-1': 0.36383246,
            'R-I-1': 1.746279178,
            'N/A': 0.00
        },
        '2.95C': {
            'total_energy': 26.99,
            'R-H-1': 10.46209293,
            'R-M-1': 0.215904013,
            'R-S-1': 15.62474914,
            'R-A-1': 0.123695901,
            'R-I-1': 0.567534062,
            'N/A': 0.00
        },
    }

    plume_data['lau_basin_abe'] = {
        '4.9C': {
            'total_energy': 19.78,
            'R-H-1': 0.382291963,
            'R-M-1': 0.274847769,
            'R-S-1': 18.33565331,
            'R-A-1': 0.518169514,
            'N/A': 0.265246444
        },
        '2.95C': {
            'total_energy': 6.49,
            'R-H-1': 0.124238398,
            'R-M-1': 0.089983413,
            'R-S-1': 6.022871649,
            'R-A-1': 0.169464288,
            'N/A': 0.083395757
        },
    }

    plume_data['lau_basin_kilo_moana'] = {
        '4.9C': {
            'total_energy': 31.74,
            'R-H-1': 0.894099271,
            'R-M-1': 0.382765378,
            'R-S-1': 28.81186728,
            'R-A-1': 0.475985621,
            'R-I-1': 1.040083982,
            'N/A': 0.137186218
        },
        '2.95C': {
            'total_energy': 10.41,
            'R-H-1': 0.290604642,
            'R-M-1': 0.125310002,
            'R-S-1': 9.465653369,
            'R-A-1': 0.156861129,
            'R-I-1': 0.330639168,
            'N/A': 0.04406072
        },
    }

    plume_data['lau_basin_mariner'] = {
        '4.9C': {
            'total_energy': 37.13,
            'R-H-1': 0.353382001,
            'R-M-1': 0.253761377,
            'R-S-1': 29.59748543,
            'R-A-1': 0.468382423,
            'R-I-1': 4.8594219,
            'R-Mn-1': 1.596469403
        },
        '2.95C': {
            'total_energy': 12.23,  # size of node, but size of node could also be proportion of reaction
            'R-H-1': 0.114585213,
            'R-M-1': 0.083011873,
            'R-S-1': 9.717413859,
            'R-A-1': 0.155712773,
            'R-I-1': 1.574375035,
            'R-Mn-1': 0.58227373
        },
    }

    plume_data_percent_yield['guaymas_basin'] = {
        '4.9C': {
            'R-H-1': 21.2,
            'R-M-1': 2.1,
            'R-S-1': 0.0,
            'R-S-2': 76.6,
            'R-I-2': 0.2,
            'N/A': 0.0
        },
        '2.95C': {
            'R-H-1': 1.3,
            'R-M-1': 92.2,
            'R-S-1': 6.4,
            'R-S-2': 0.0,
            'R-I-2': 0.0,
            'N/A': 0.0
        },
    }

    plume_data_percent_yield['cayman_von_damm'] = {
        '4.9C': {
            'R-H-1': 50.6,
            'R-M-1': 25.9,
            'R-S-1': 16.3,
            'R-S-2': 0.8,
            'R-A-1': 5.9,
            'N/A': 0.6
        },
        '2.95C': {
            'R-H-1': 51.1,
            'R-M-1': 26.4,
            'R-S-1': 21.8,
            'R-S-2': 0.0,
            'R-A-1': 0.0,
            'N/A': 0.7
        },
    }

    plume_data_percent_yield['cayman_picard'] = {
        '4.9C': {
            'R-H-1': 39.1,
            'R-M-1': 0.8,
            'R-S-1': 57.6,
            'R-A-1': 0.4,
            'R-I-1': 2.1,
            'N/A': 0.0
        },
        '2.95C': {
            'R-H-1': 38.8,
            'R-M-1': 0.8,
            'R-S-1': 57.9,
            'R-A-1': 0.5,
            'R-I-1': 2.1,
            'N/A': 0.0
        },
    }

    plume_data_percent_yield['lau_basin_abe'] = {
        '4.9C': {
            'R-H-1': 1.9,
            'R-M-1': 1.4,
            'R-S-1': 92.7,
            'R-A-1': 2.6,
            'N/A': 1.3
        },
        '2.95C': {
            'R-H-1': 1.9,
            'R-M-1': 1.4,
            'R-S-1': 92.8,
            'R-A-1': 2.6,
            'N/A': 1.3
        },
    }

    plume_data_percent_yield['lau_basin_kilo_moana'] = {
        '4.9C': {
            'R-H-1': 2.8,
            'R-M-1': 1.2,
            'R-S-1': 90.8,
            'R-A-1': 1.5,
            'R-I-1': 3.3,
            'N/A': 0.4
        },
        '2.95C': {
            'R-H-1': 2.8,
            'R-M-1': 1.2,
            'R-S-1': 90.9,
            'R-A-1': 1.5,
            'R-I-1': 3.2,
            'N/A': 0.4
        },
    }

    plume_data_percent_yield['lau_basin_mariner'] = {
        '4.9C': {
            'R-H-1': 1.0,
            'R-M-1': 0.7,
            'R-S-1': 79.7,
            'R-A-1': 1.3,
            'R-I-1': 13.1,
            'R-Mn-1': 4.3
        },
        '2.95C': {
            'R-H-1': 0.9,
            'R-M-1': 0.7,
            'R-S-1': 79.5,
            'R-A-1': 1.3,
            'R-I-1': 12.9,
            'R-Mn-1': 4.8
        },
    }

    # -------------- Make into networks ----------------

    # turn these into a directed bipartite network with edge attributes
    # bipartite graph
    # reaction id is the reaction node
    # two graphs per plume, one for 4.9C and the other for 2.95C

    plume_networks = {}

    for plume in plumes:

        G = nx.DiGraph()

        for temp in temps:

            plume_reaction_energies = plume_data_percent_yield[plume][temp]

            for reaction_id in plume_reaction_energies.keys():

                if reaction_id == 'total_energy' or reaction_id == 'N/A':
                    continue

                reaction_data_part = reaction_data[reaction_id]
                electrons = reaction_data_part['electrons']
                dg_1 = reaction_data_part['dg_1']
                dg_25 = reaction_data_part['dg_25']
                dg_100 = reaction_data_part['dg_100']

                # bipartite graph with reactions as one node color
                G.add_node(reaction_id, bipartite=1, electrons=electrons, dg_1=dg_1, dg_25=dg_25, dg_100=dg_100, energy=plume_reaction_energies[reaction_id])

                # in edges and nodes for reaction node
                for in_node in reaction_data_part['ins']:

                    # look for amount edge attribute
                    if re.search('^[\d\.]+', in_node):
                        amount = re.findall('^[\d\.]+', in_node)[0]
                        chem = in_node[len(amount):]
                        amount = float(amount)
                    else:
                        amount = 1.0
                        chem = in_node

                    G.add_node(chem, bipartite=0)
                    G.add_edge(chem, reaction_id, amount=amount)

                # out edges and nodes for reaction node
                for out_node in reaction_data_part['outs']:

                    # look for amount edge attribute
                    if re.search('^[\d\.]+', out_node):
                        amount = re.findall('^[\d\.]+', out_node)[0]
                        chem = out_node[len(amount):]
                        amount = float(amount)
                    else:
                        amount = 1.0
                        chem = out_node

                    G.add_node(chem, bipartite=0)
                    G.add_edge(reaction_id, chem, amount=amount)

            plume_networks[plume+'_'+temp] = G

    # apply BDM measure to each graph
    # perturb each reaction node to see the change in complexity
    # save that value as a node attribute

    bdm = BDM(ndim=2)

    for plume_network in plume_networks.keys():

        # get bdm for graph
        G = plume_networks[plume_network]
        bottom_nodes, top_nodes = bipartite.sets(G)

        adj = nx.adj_matrix(G)
        adj = scipy.sparse.csr_matrix.toarray(adj)
        bdm_before = bdm.bdm(adj)

        for reaction_node in bottom_nodes:

            # get bdm for graph with removed node
            removed_graph = copy.deepcopy(G)
            removed_graph.remove_node(reaction_node)
            adj = nx.adj_matrix(removed_graph)
            adj = scipy.sparse.csr_matrix.toarray(adj)
            bdm_after = bdm.bdm(adj)

            # add difference in node attribute
            G.nodes[reaction_node]['dC'] = bdm_before - bdm_after
            G.nodes[reaction_node]['%dC'] = (bdm_before - bdm_after)/bdm_before

            # add network structure measures here
            G.nodes[reaction_node]['degree'] = list(G.degree([reaction_node]))[0][1]

    # -------------- Plot results ----------------

    # scatterplot of complexity vs. one of (degree, electrons, dg_1, dg_25, dg_100, energy)
    # get node attributes into a df

    node_data = []

    for network in plume_networks.keys():

        node_dict = plume_networks[network].nodes._nodes
        node_dict2 = {}
        for node in node_dict.keys():
            nd = node_dict[node]
            if nd['bipartite'] == 1:
                nd['Reaction ID'] = node
                node_dict2[node] = nd

        # make into a df to concat together with others
        node_df = pd.DataFrame.from_dict(node_dict2, orient='index')
        plume = re.findall('[^\d]+', network)[0][:-1]
        temp = network.split('_')[-1]
        node_df['Plume'] = plume
        node_df['Sample Temp'] = temp
        node_data.append(node_df)

    node_data = pd.concat(node_data)

    # scatterplots
    sns.set(rc={'figure.figsize': (10, 7)})
    sns.set(font_scale=1.5)
    scatter = sns.scatterplot(data=node_data, x="dC", y="energy", hue="Plume", style="Reaction ID", s=100, linewidth=0)
    #scatter = sns.scatterplot(data=node_data, x="degree", y="energy", hue="Network")
    #scatter.legend(fontsize=6)
    plt.axvline(0, color='grey', linestyle='--', linewidth=1)

    # add trendline for whole thing

    X_plot = node_data["dC"]
    Y_plot = node_data["energy"]
    m, b = np.polyfit(X_plot, Y_plot, 1)
    plt.plot(X_plot, m*X_plot + b, color='r')

    plt.ylabel('% energy yield')

    sns.despine()
    scatter.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=12)
    plt.tight_layout()
    plt.savefig('dc_vs_percent_energy.pdf')
    #plt.savefig('degree_vs_energy.pdf')
    plt.show()
    plt.clf()

    X_plot = node_data["dC"]
    Y_plot = node_data["energy"]
    X2 = sm.add_constant(X_plot)
    est = sm.OLS(Y_plot, X2)
    est2 = est.fit()
    print("summary()\n", est2.summary())
    print("pvalues\n", est2.pvalues)
    print("f_pvalues\n", est2.f_pvalue)
    print("tvalues\n", est2.tvalues)
    print("rsquared\n", est2.rsquared)
    print("rsquared_adj\n", est2.rsquared_adj)

    # visualize networks

    for plume_network in plume_networks.keys():

        G = plume_networks[plume_network]
        bottom_nodes, top_nodes = bipartite.sets(G)

        color_dict = {0: 'silver', 1: 'y'}
        color_list = [color_dict[i[1]] for i in G.nodes.data('bipartite')]

        # Draw bipartite graph
        #pos = dict()
        color = []
        #pos.update((n, (1, i)) for i, n in enumerate(bottom_nodes))  # put nodes from X at x=1
        #pos.update((n, (2, i)) for i, n in enumerate(top_nodes))  # put nodes from Y at x=2
        pos = nx.spring_layout(G)

        # node size scale with reaction energy (only reaction nodes)
        node_attributes = dict(G.nodes.data())
        nodes = list(G.nodes())
        node_sizes = []
        scale = 10
        for node in nodes:
            try:
                node_sizes.append(scale*node_attributes[node]['energy'])
            except:
                node_sizes.append(10)

        # scaled nodes
        #nx.draw(G, pos=pos, with_labels=True, node_color=color_list, font_size=8, node_size=node_sizes)
        nx.draw(G, pos=pos, with_labels=True, node_color=color_list, font_size=8, node_size=node_sizes)
        plt.show()
        plt.savefig(plume_network + '.pdf')
        plt.clf()

