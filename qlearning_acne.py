import pandas as pd
import numpy as np
import argparse
import time
from env import GeneEnv
from agent import Agent
import matplotlib.pyplot as plt

# Argument parser to choose dataset
parser = argparse.ArgumentParser(description="Q-Learning Model for Acne Gene Selection")
parser.add_argument('--dataset', type=str, required=True, help="Dataset name (e.g., GSE6475, GSE53795, GSE108110)")
args = parser.parse_args()

# File paths based on dataset argument
dataset = args.dataset
data_file = f'./data/data_degenes_with_acne_sample_{dataset}.csv'
ggc_file = f'./data/WGCNA_Matrix_{dataset}_Final.csv'
top_gene_file = f'./data/gene_sign_list_{dataset}.csv'

# Load data files
df = pd.read_csv(data_file, index_col=0)
ggc_df = pd.read_csv(ggc_file, index_col=0)
top_gene_list = pd.read_csv(top_gene_file)

# Setup parameters
samples = [np.array(df.iloc[i].values.tolist()) for i in range(6)]
sample_obj = {f"sample{i+1}": sample for i, sample in enumerate(samples)}
# temp_gene_id = ggc_df.index
gene_id = ggc_df.index.tolist()
top_gene_id_list = top_gene_list.iloc[:,0].tolist()
count_row = len(ggc_df)
num_states = 2**len(top_gene_list)
num_actions = count_row
terminated_index = count_row

num_episodes = 2000
record = 0
total_score = 0
final = []
selected_gene_counter = []

plot_connectivity = []
plot_mean_scores = []
plot_rewards = []

st = time.time()
timestamp = time.strftime("%Y%m%d-%H%M%S")

# Training loop
for key, value in sample_obj.items():
    print(f"Processing {key} for dataset {dataset}...")
    env = GeneEnv(ggc_df, top_gene_id_list, count_row, gene_id, value)
    agent = Agent(num_states, num_actions, num_episodes, terminated_index)
    top_gene_list_index = env.get_index_top_genes()
    
    for episodes in range(num_episodes):
        terminated = False
        log_action = [terminated_index]
        cul_reward = 0
        
        while not terminated:
            state_old = agent.get_state(env, terminated_index, top_gene_list_index)
            state_old_index = agent.get_state_index(state_old)
            action = agent.get_action(state_old_index)
            
            if action not in log_action:
                terminated, selection, gene_expressed = env.play_step(action, terminated_index)
                log_action.append(action)
                
                selected_gene = agent.get_selected_gene(action)
                reward = agent.cal_reward(log_action, env, selected_gene, top_gene_list_index, terminated, gene_expressed)
                cul_reward += reward
                
                state_new = agent.get_state(env, terminated_index, top_gene_list_index)
                state_new_index = agent.get_state_index(state_new)
                agent.update_Qtable(state_old_index, action, reward, state_new_index, terminated)
            else:
                action = terminated_index
                terminated, selection, gene_expressed = env.play_step(action, terminated_index)
                log_action.append(action)
                reward = agent.cal_reward(log_action, env, selected_gene, top_gene_list_index, terminated, gene_expressed)
                cul_reward += reward
                state_new = agent.get_state(env, terminated_index, top_gene_list_index)
                state_new_index = agent.get_state_index(state_new)
                agent.update_Qtable(state_old_index, action, reward, state_new_index, terminated)
        
        plot_rewards.append(cul_reward)
        connectivity = env.cal_max_connectivity(selection)
        if connectivity > record:
            record = connectivity
                
        plot_connectivity.append(connectivity)
        total_score += connectivity
        plot_mean_scores.append(total_score / (episodes + 1))
        
        print(f"Episode {episodes}: Cumulative Reward: {cul_reward}, Connectivity: {connectivity}, Record: {record}")
        selected_gene_counter.append(log_action)
        env.reset()

# Plotting results
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(plot_rewards, label='Cumulative Reward')
plt.xlabel('Episodes')
plt.ylabel('Reward')
plt.title('Cumulative Reward Over Episodes')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(plot_connectivity, label='Connectivity', color='orange')
plt.plot(plot_mean_scores, label='Mean Score', color='green')
plt.xlabel('Episodes')
plt.ylabel('Scores')
plt.title('Connectivity and Mean Scores Over Episodes')
plt.legend()

plt.tight_layout()
plt.savefig(f'training_{dataset}_results_plot_{timestamp}.png')

# Exporting results
max_length = max(len(row) for row in selected_gene_counter)
padded_selected_gene_counter = [row + [None] * (max_length - len(row)) for row in selected_gene_counter]
selected_gene_df = pd.DataFrame(padded_selected_gene_counter, columns=[f"Gene{idx+1}" for idx in range(max_length)])
selected_gene_df.to_csv(f'final_selected_gene_count_{dataset}_{num_episodes}.csv', index=False)

q_table = agent.model.get_Q_table()
q_table_df = pd.DataFrame(q_table)
q_table_df.to_csv(f'final_q_table_{dataset}_{num_episodes}.csv')

et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
