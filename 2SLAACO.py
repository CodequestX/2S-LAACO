import operator
import networkx as nx
import csv
import random
import matplotlib.pyplot as plt
import copy
import time
import collections
import itertools
import math
import networkx.algorithms.traversal.breadth_first_search as bfs
import numpy as np
import time

start_time = time.time()

#function for computing intra-group connectedness
def get_intra_group_connectedness(group, gr):
    d = 0
    sg = gr.subgraph(group)
    if (len(group) - 1) > 0:
        d = (len(sg.edges()))/(len(group) * (len(group) - 1))
    else:
        d = 1
    
    return d

#function for computing inter-group connectedness
def get_inter_group_connectedness(group, gr):
    d = 0
    sg = gr.subgraph(group)
    if len(gr.edges(group)) > 0:
        d = (len(gr.edges(group)) - len(sg.edges()))/len(gr.edges(group))
    else:
        d = 1
    
    return d

#function for computing a node's score using the MGC-c heuristic 
def get_contribution(theta, node, gr, group_mapping, groups, intra_group_conn, inter_group_conn):
    contr = 0
    nbrs = gr.neighbors(node)
    grp_contr = {}
    grp_contr_num = {}

    grp_coverage = []
    if node in group_mapping.keys():
        grp_coverage = group_mapping[node]

    
    for xg in range(len(grp_coverage)):
        grp_size =  len(groups[xg])
        grp_contr[xg] = 1/(theta*grp_size)
        
    for n in nbrs:
        for gx in grp_coverage:
            if n in group_mapping.keys():
                if gx in group_mapping[n]:
                    grp_size =  len(groups[int(gx)])
                    if gx in grp_contr.keys():
                        grp_contr[gx] = grp_contr[gx] + 1/(theta*grp_size)
                    else:
                        grp_contr[gx] =  1/(theta*grp_size)
    for g in grp_contr:
        contr = contr + grp_contr[g] * (intra_group_conn[int(g)] + inter_group_conn[int(g)])
    return contr

#function for calculating node activation probabilities 
def calc_act_prob(lg, node):
    nodes = lg.nodes(data=True)
    parent = 0
    act_prev = nodes[node]['act']
    xpp = 1
    for pred in lg.predecessors(str(node)):
        parent = nodes[pred]['act']
        p_uv = 1 - pow((1 - .01), math.fabs(lg.get_edge_data(str(pred), str(node))["weight"]))
        xpp = xpp * (1 - math.fabs(p_uv)*parent)
            
    act_t = (1 - xpp) * (1 - act_prev)
    act = act_prev + act_t
    
    return act_t, act

#function for calculating group influence spread
def calc_IS(gr, seedset, groups, theta):
    #lg = nx.DiGraph()
    #all_node_list_by_level = []
    infl_spread_t ={}
    for i in gr.nodes(): 
        infl_spread_t[i] = 0
    infl_spread = [0] * len(groups)
    t = 0
    active_groups = 0
    
    for v in gr.nodes():
        gr.nodes[str(v)]['act'] = 0
        gr.nodes[str(v)]['act_t'] = 0
        if v in seedset:
            gr.nodes[str(v)]['act'] = 1
            gr.nodes[str(v)]['act_t'] = 1
    
            
    while(t < 9):
        for v in gr.nodes():
            act_t, act = calc_act_prob(gr, v)
            infl_spread_t[v] = infl_spread_t[v] + act_t
            gr.nodes[str(v)]['act'] = act
            gr.nodes[str(v)]['act_t'] = act_t

        t = t + 1

    gc = 0
    for g in groups:
        for v in g:
            infl_spread[gc] = infl_spread[gc] + gr.nodes[str(v)]['act']
        if infl_spread[gc] > theta * len(g):
            active_groups = active_groups + 1
        gc = gc + 1
        
    return active_groups

#function to construct a candidate seed set for each ant
def construct_candidate_seed_set(k, candidate_pool, ant):
    selected_nodes = set()
    for _ in range(k):
        node = select_node(selected_nodes, candidate_pool, ant)
        selected_nodes.add(node)
    return selected_nodes

#function to select the next node for an ant based on pheromone levels and MGC-c heuristic
def select_node(selected_nodes, candidate_pool, ant):
    probabilities = []
    for node in candidate_pool:
        if node not in selected_nodes:
            probability = (pheromones[node] ** alpha[ant])  * (node_contr_top[node] ** beta[ant])
            probabilities.append((node, probability))
    
    #normalize the probabilities
    total_probability = sum(prob for _, prob in probabilities)
    probabilities = [(node, prob / total_probability) for node, prob in probabilities]

    selected_node = random.choices([node for node, _ in probabilities], [prob for _, prob in probabilities])[0]
    
    return selected_node

#function to update pheromones after each iteration
def update_pheromones(all_candidate_seed_sets, all_influence_spreads, candidate_pool):
    global pheromones
    #evaporate pheromones
    for node in candidate_pool:
        pheromones[node] *= (1 - evaporation_rate)
    #identify the candidate seed sets with high group influence spread
    
    print('all_candidate_seed_sets')
    print(all_candidate_seed_sets)
    print('all_influence_spreads')
    print(all_influence_spreads)
    
    candidate_seed_sets_inf_list = list(zip(all_candidate_seed_sets, all_influence_spreads)) 
    candidate_seed_sets_inf_list.sort(key=lambda x: x[1])
    #candidate_seed_sets_inf_sorted = sorted(candidate_seed_sets_inf_dict.items(), key=operator.itemgetter(1), reverse = True)
    high_inf_candidate_seed_sets = candidate_seed_sets_inf_list[0:math.ceil(.25*num_ants)]
    high_inf_candidate_seed_set_list = [t[0] for t in high_inf_candidate_seed_sets] 
    
    #deposit pheromone based on group influence spread
    for ant_candidate_seed_set, influence in zip(all_candidate_seed_sets, all_influence_spreads):
        #pheromone_deposit = influence
        for seed_node in ant_candidate_seed_set:
            if ant_candidate_seed_set in high_inf_candidate_seed_set_list:
                pheromones[seed_node] += influence
                

#function to select an action based on current probabilities
def select_action(num_actions_LA, action_probabilities_LA):
    return random.choices(range(num_actions_LA), action_probabilities_LA)[0]


#function to generate feedback
def get_feedback(param_name, ant, num_ants, selected_action, all_influence_spreads_prev, all_influence_spreads):
    if len(all_influence_spreads_prev) > 0:
        if param_name == 'evaporation_rate':
            rand_ant1 = random.randint(0,num_ants)
            rand_ant2 = random.choice([i for i in range(0,num_ants) if i not in [rand_ant1]])
            
            avg_influence_spread_prev = sum(all_influence_spreads_prev) / num_ants
            avg_influence_spread = sum(all_influence_spreads) / num_ants
            
            if avg_influence_spread <= avg_influence_spread_prev or all_influence_spreads[rand_ant1] == all_influence_spreads[rand_ant2]:
                feedback = 1
            else:
                feedback = 0
        else:
            rand_ant = random.choice([i for i in range(0,num_ants) if i not in [ant]])
            if all_influence_spreads[ant] <= all_influence_spreads_prev[ant] or all_influence_spreads[ant] == all_influence_spreads[rand_ant]:
                feedback = 1
            else:
                feedback = 0
    else:
        feedback = 0
    return feedback

#function to update action probabilities based on feedback
def update_action_probabilities_LA(selected_action, action_probabilities_LA, num_actions_LA, feedback, reward_factor_LA, penalty_factor_LA):
    for action in range(num_actions_LA):
        if feedback == 0:  #favourable feedback
            if action == selected_action:
                #increase the probability of the successful action
                action_probabilities_LA[action] = action_probabilities_LA[action] + reward_factor_LA * (1 - action_probabilities_LA[action])
            else:
                #decrease the probabilities of the other actions
                action_probabilities_LA[action] = (1 - penalty_factor_LA) * action_probabilities_LA[action]
        else:  #unfavourable feedback
            if action == selected_action:
                #decrease the probability of the unsuccessful action
                action_probabilities_LA[action] = (1 - reward_factor_LA) * action_probabilities_LA[action]
                #increase the probability of the other actions
                action_probabilities_LA[action] = (penalty_factor_LA / (num_actions_LA - 1)) + (1 - penalty_factor_LA) * action_probabilities_LA[action]
                
    #ensure that the action probabilities are between 0 and 1
    action_probabilities_LA = [max(0.0, min(prob, 1.0)) for prob in action_probabilities_LA]
    return action_probabilities_LA


#function to update parameter based on selected action
def update_parameter(param, selected_action, delta):
    if selected_action == 0:
        param = param + delta
    if selected_action == 1:
        param = param - delta
    return param


#function for the LA operations
def LA_operations(param_name, ant, num_ants, curr_selected_action, param, action_probabilities_LA, delta, num_actions_LA, reward_factor_LA, penalty_factor_LA, all_influence_spreads_prev, all_influence_spreads):
    #update the parameter based on current selected action
    param = update_parameter(param, curr_selected_action, delta) 
    #get feedback from the environment for the current selected action
    feedback = get_feedback(param_name, ant, num_ants, curr_selected_action, all_influence_spreads_prev, all_influence_spreads)

    #update the action probabilities based on feedback
    action_probabilities_LA = update_action_probabilities_LA(curr_selected_action, action_probabilities_LA, num_actions_LA, feedback, reward_factor_LA, penalty_factor_LA)
    
        
    #select the next action based on the updated action probabilities
    next_selected_action = select_action(num_actions_LA, action_probabilities_LA)
    
    return next_selected_action, param


#function for GIM ant colony optimization with adaptive parameter setting using LA
def ant_colony_optimization_LA(k, num_ants, num_iterations, candidate_pool, action_probabilities_LA_evaporation_rate, action_probabilities_LA_alpha, action_probabilities_LA_beta, evaporation_rate, alpha, beta, delta, num_actions_LA):
    best_candidate_seed_set = None
    best_influence = 0
    
    #initialize LA actions
    selected_action_evaporation_rate = 2
    selected_action_alpha = [-1]*num_ants
    selected_action_beta = [-1]*num_ants
    for ant in range(num_ants):
        selected_action_alpha[ant] = 2
        selected_action_beta[ant] = 2
        
    all_influence_spreads = []
    all_influence_spreads_prev = []            
    for iteration in range(num_iterations):
        print(iteration)
        all_influence_spreads_prev = all_influence_spreads
        all_candidate_seed_sets = []
        all_influence_spreads = []
         
        
        #each ant constructs a solution (selects nodes)
        for ant in range(num_ants):
            ant_candidate_seed_set = construct_candidate_seed_set(k, candidate_pool, ant)
            influence = calc_IS(gr, ant_candidate_seed_set, groups, theta)
            
            all_candidate_seed_sets.append(ant_candidate_seed_set)
            all_influence_spreads.append(influence)
            
            if influence > best_influence:
                best_influence = influence
                best_candidate_seed_set = ant_candidate_seed_set
                
        
        #update pheromone levels after all ants have selected their nodes
        update_pheromones(all_candidate_seed_sets, all_influence_spreads, candidate_pool)
        
        #adaptive parameter setting using LA
        evaporation_rate = LA_operations('evaporation_rate', -1, num_ants, selected_action_evaporation_rate, evaporation_rate, action_probabilities_LA_evaporation_rate, delta, num_actions_LA, reward_factor_LA, penalty_factor_LA, all_influence_spreads_prev, all_influence_spreads)
        for ant in range(num_ants):
            selected_action_alpha[ant], alpha[ant] = LA_operations('alpha', ant, num_ants, selected_action_alpha[ant], alpha[ant], action_probabilities_LA_alpha[ant], delta, num_actions_LA, reward_factor_LA, penalty_factor_LA, all_influence_spreads_prev, all_influence_spreads)
            selected_action_beta[ant], beta[ant] = LA_operations('beta', ant, num_ants, selected_action_beta[ant], beta[ant], action_probabilities_LA_beta[ant], delta, num_actions_LA, reward_factor_LA, penalty_factor_LA, all_influence_spreads_prev, all_influence_spreads)
            
        
    return best_candidate_seed_set, best_influence
    

#initialization    
gr = nx.DiGraph() #graph
theta = 0.5 #activation threshold
num_iterations = 150 #maximum number of iterations
num_ants = 50 #number of ants in the ant colony
alpha = [0.2] * num_ants #list of initial value of alpha for each ant
beta = [0.8] * num_ants #list of initial value of beta for each ant
evaporation_rate = 0.5 #initial pheromone evaporation rate
div = 0.3 #diversity factor
avg_grp_size = 0 #initialize average group size
k = 10 #length of seed user set

num_actions_LA = 3 #number of actions for LA
action_probabilities_LA_evaporation_rate = [1.0 / num_actions_LA] * num_actions_LA #LA action probabilities for evaporation rate (all actions start with an equal probability)
action_probabilities_LA_alpha = [[1.0 / num_actions_LA] * num_actions_LA]*num_ants #LA action probabilities for alpha for each ant (all actions start with an equal probability)
action_probabilities_LA_beta = [[1.0 / num_actions_LA] * num_actions_LA]*num_ants #LA action probabilities for beta for each ant (all actions start with an equal probability)
reward_factor_LA = 0.01 #reward factor for adjusting LA action probabilities
penalty_factor_LA =0.01 #penalty factor for adjusting LA action probabilities
delta = 0.1 #step value
#build the graph from csv
with open(r"dataset.csv") as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        gr.add_nodes_from([row[0], row[1]], shapley = 0)
        gr.add_edge(row[0], row[1], weight = float(row[2]))

#read the groups
group_mapping = {}

with open(r"group.csv") as csvDataFile:
        csvReader = csv.reader(csvDataFile)
        count = 0
        for row in csvReader:
            for i in range(len(row)):
                if row[i] in group_mapping.keys():
                    group_mapping[row[i]].append(count)
                else:
                    group_mapping[row[i]] = [count]
            count = count + 1
            
groups = []
with open(r"group.csv") as csvDataFile:
        csvReader = csv.reader(csvDataFile)
        for row in csvReader:
                groups.append(row)
                
#calculate intra group connectedness for all groups
intra_group_conn = []
for g in range(len(groups)):
    intra_group_conn.append(get_intra_group_connectedness(groups[int(g)], gr))
    
#calculate intra group connectedness for all groups
inter_group_conn = []
for g in range(len(groups)):
    inter_group_conn.append(get_inter_group_connectedness(groups[int(g)], gr))

#calculate the average size of groups
sum_grp_size = 0
for c in range(len(groups)):
    sum_grp_size = sum_grp_size + len(groups[c])
avg_grp_size = sum_grp_size/len(groups)
    
nodeList = gr.nodes()


#-------------------------------------begin:stage 1----------------------------------------
#compute pool size
pool_size = math.ceil(k + (len(nodeList)-k) * pow(((div*k*theta*avg_grp_size) / len(nodeList)), (1 - div)))

#construct the candidate pool
node_contr_sorted_dict = {}
nodeList = list(gr.nodes())

node_contr = {}
for node in nodeList:
    node_contr[node] = get_contribution(theta, node, gr, group_mapping, groups, intra_group_conn, inter_group_conn)
node_contr_sorted = {}
node_contr_sorted = sorted(node_contr.items(), key=operator.itemgetter(1), reverse = True)

node_contr_top = dict(node_contr_sorted[0:pool_size])
candidate_pool = list(node_contr_top.keys())
#-------------------------------------end:stage 1----------------------------------------

#-------------------------------------begin:stage 2----------------------------------------
#Initialize pheromone dictionary with each node's pheromone value set to 1
pheromones = {node: 1.0 for node in candidate_pool}
best_candidate_seed_set, best_influence = ant_colony_optimization_LA(k, num_ants, num_iterations, candidate_pool, action_probabilities_LA_evaporation_rate, action_probabilities_LA_alpha, action_probabilities_LA_beta, evaporation_rate, alpha, beta, delta, num_actions_LA)
#-------------------------------------end:stage 2----------------------------------------

print(best_candidate_seed_set)
print(best_influence)
print("--- %s seconds ---" % (time.time() - start_time))
