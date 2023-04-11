import pandas as pd
import numpy as np
from env import GeneEnv
from agent import Agent
from helper import plot
import matplotlib.pyplot as plt
import time

st = time.time()
df = pd.read_csv('./data/data_degenes_with_ance_sample_gse53795.csv', index_col=0)
ggc_df = pd.read_csv('./data/WGCNA_Matrix_GSE53795_Final.csv', index_col=0)
top_gene_list = pd.read_csv('./data/gene_sign_list_GSE53795.csv')

sample1 = np.array(df.iloc[0].values.tolist())
sample2 = np.array(df.iloc[1].values.tolist())
sample3 = np.array(df.iloc[2].values.tolist())
sample4 = np.array(df.iloc[3].values.tolist())
sample5 = np.array(df.iloc[4].values.tolist())
sample6 = np.array(df.iloc[5].values.tolist())
sample7 = np.array(df.iloc[6].values.tolist())
sample8 = np.array(df.iloc[7].values.tolist())
sample9 = np.array(df.iloc[8].values.tolist())
sample10 = np.array(df.iloc[9].values.tolist())
sample11 = np.array(df.iloc[10].values.tolist())
sample12 = np.array(df.iloc[11].values.tolist())

sample_obj = {
    "sample1" : sample1,
    "sample2" : sample2,
    "sample3" : sample3,
    "sample4" : sample4,
    "sample5" : sample5,
    "sample6" : sample6,
    "sample7" : sample7,
    "sample8" : sample8,
    "sample9" : sample9,
    "sample10" : sample10,
    "sample11" : sample11,
    "sample12" : sample12
}

temp_gene_id = ggc_df.index
gene_id = temp_gene_id.tolist()


temp_top_gene_id = top_gene_list.iloc[:,0]
top_gene_id_list = temp_top_gene_id.tolist()
print("top_gene_list: ", top_gene_id_list)
print("gene id: ",gene_id)
print("ggc_df: ",ggc_df)

count_row = len(ggc_df)
print(count_row)

size = len(top_gene_list)
    
num_states = 2**(10)
print(num_states)
state_space = []
num_actions = count_row

terminated_index = count_row # int for loop break, index that out of range
log_action = [terminated_index] # index that make the read from array out of range, to represent the terminated
terminated = False

n = 0
num_iteration = 0
num_episodes = 200
record = 0
total_score = 0
plot_connectivity = []
plot_mean_scores = []
plot_rewards = []

cul_reward = 0
reward=0
selected_gene = 0
selected_gene_list = []
selected_gene_lists = []
gene_markers = []
final = []
actions_space = num_actions+1
q_table = np.zeros((num_states, actions_space))

for key, value in sample_obj.items() :
    print("sample: ",key, value)
    env = GeneEnv(ggc_df, top_gene_id_list, count_row, gene_id, value)
    agent = Agent(num_states, num_actions, num_episodes, terminated_index)
    top_gene_list_index = env.get_index_top_genes()
# Start training
    for episodes in range(num_episodes):
        num_iteration = num_iteration + 1
        print("loop ", episodes)
        print("check reward: ", reward)
        cul_reward = 0
        while (not(terminated)):
            state_old = agent.get_state(env, terminated_index, top_gene_list_index) # need to convert to index
            state_old_index = agent.get_state_index(state_old)
            print("old state is: ", state_old)
            print("old state index is: ", state_old_index)

            action = agent.get_action(state_old_index) # use state index to get action
            print("action : ", action)

            ## check is the action have been taken?
            ## if not, then is a valid action; if yes, then terminated
            # action = selected gene
            if (action not in log_action):
                terminated, selection, gene_expressed = env.play_step(action, terminated_index)
                log_action.append(action)

                selected_gene = agent.get_selected_gene(action)
                selected_gene_list.append(selected_gene)
                selected_gene_lists.append(selected_gene)
                print("len(selected_gene_list) :", len(selected_gene_list))

                reward = agent.cal_reward(log_action, env, selected_gene, top_gene_list_index, terminated, gene_expressed)
                cul_reward += reward
                state_new = agent.get_state(env, terminated_index, top_gene_list_index)
                state_new_index = agent.get_state_index(state_new)
                agent.update_Qtable(state_old_index, action, reward, state_new_index, terminated)       
            else:
                print("action has been repeatedlly selected and considered action index as terminated_index")
                action = terminated_index
                terminated, selection, gene_expressed = env.play_step(action, terminated_index)
                log_action.append(action)

                selected_gene = agent.get_selected_gene(action)
                selected_gene_list.append(selected_gene)
                selected_gene_lists.append(selected_gene)

                reward = agent.cal_reward(log_action, env, selected_gene, top_gene_list_index, terminated, gene_expressed)
                cul_reward += reward
                state_new = agent.get_state(env, terminated_index, top_gene_list_index)
                state_new_index = agent.get_state_index(state_new)
                log_action= [terminated_index]

                agent.update_Qtable(state_old_index, action, reward, state_new_index, terminated)
                plot_mean_scores.append(reward)
                reward = 0

        plot_rewards.append(cul_reward)
        connectivity = env.cal_max_connectivity(selection)
        if connectivity > record:
            record = connectivity                
        print("culmulative reward :", cul_reward)
        print('Number of episodes', episodes, 'Connectivity', connectivity, 'Record:', record)

        plot_connectivity.append(cul_reward)
        total_score += connectivity
        mean_score = total_score / (episodes+1)
        
        plot(plot_connectivity, plot_mean_scores, num_iteration)
        print("cul reward:", cul_reward)
        terminated = False
        
        if (num_iteration >= 2390):
            final = selected_gene_list
            q_table = agent.model.get_Q_table()
        selected_gene_list = []
        env.reset()

gene_markers_count = []
print("final :",final)

et = time.time()
elapsed_time = et - st

for item in (selected_gene_lists):
    if ((item != terminated_index) and (gene_id[item] not in gene_markers)):
        count = selected_gene_lists.count(item)
        gene_markers.append(gene_id[item])
        gene_markers_count.append(count)

selected_gene_count = zip(gene_markers, gene_markers_count)
selected_gene_count_df = pd.DataFrame(selected_gene_count)
selected_gene_count_df.to_csv('final_v2_cross_sample_selected_gene_count_stateGSE53795_2400.csv')
print(plot_rewards)
print('Execution time:', elapsed_time, 'seconds')
q_table_df = pd.DataFrame(q_table)
q_table_df.to_csv('final_v2_cross_sample_q-table_2400_gse53795.csv')
