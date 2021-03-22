import csv
import re
import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

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

    plume_data = {}

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
            'total_energy': 12.23,
            'R-H-1': 0.114585213,
            'R-M-1': 0.083011873,
            'R-S-1': 9.717413859,
            'R-A-1': 0.155712773,
            'R-I-1': 1.574375035,
            'R-Mn-1': 0.58227373
        },
    }

    # turn these into a directed bipartite network with edge attributes
    # bipartite graph
    # reaction id is the reaction node
    # two graphs per plume, one for 4.9C and the other for 2.95C

    plume_networks = {}

    for plume in plumes:

        G = nx.DiGraph()

        for temp in temps:

            plume_reaction_energies = plume_data[plume][temp]

            for reaction_id in plume_reaction_energies.keys():

                if reaction_id == 'total_energy' or reaction_id == 'N/A':
                    continue

                reaction_data_part = reaction_data[reaction_id]
                electrons = reaction_data_part['electrons']
                dg_1 = reaction_data_part['dg_1']
                dg_25 = reaction_data_part['dg_25']

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

    # plot the networks
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
        scale = 700
        for node in nodes:
            try:
                node_sizes.append(scale*node_attributes[node]['energy'])
            except:
                node_sizes.append(500)

        # scaled nodes
        #nx.draw(G, pos=pos, with_labels=True, node_color=color_list, font_size=8, node_size=node_sizes)
        nx.draw(G, pos=pos, with_labels=True, node_color=color_list, font_size=8)
        #plt.show()
        plt.savefig(plume_network + '.pdf')
        plt.clf()

