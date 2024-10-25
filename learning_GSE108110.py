import pandas as pd
import numpy as np
from env import GeneEnv
from agent import Agent
import matplotlib.pyplot as plt
import time

st = time.time()
timestamp = time.strftime("%Y%m%d-%H%M%S")

df = pd.read_csv('./data/data_degenes_with_ance_sample_gse108110.csv', index_col=0)
ggc_df = pd.read_csv('./data/WGCNA_Matrix_GSE108110_Final.csv', index_col=0)
top_gene_list = pd.read_csv('./data/gene_sign_list_GSE108110.csv')

# sample1 = np.array(df.iloc[0].values.tolist())
# sample2 = np.array(df.iloc[1].values.tolist())
# sample3 = np.array(df.iloc[2].values.tolist())
# sample4 = np.array(df.iloc[3].values.tolist())
# sample5 = np.array(df.iloc[4].values.tolist())
# sample6 = np.array(df.iloc[5].values.tolist())
# sample7 = np.array(df.iloc[6].values.tolist())
# sample8 = np.array(df.iloc[7].values.tolist())
# sample9 = np.array(df.iloc[8].values.tolist())
samples = [np.array(df.iloc[i].values.tolist()) for i in range(9)]

# sample_obj = {
#     "sample1" : sample1,
#     "sample2" : sample2,
#     "sample3" : sample3,
#     "sample4" : sample4,
#     "sample5" : sample5,
#     "sample6" : sample6,
#     "sample7" : sample7,
#     "sample8" : sample8,
#     "sample9" : sample9
# }
sample_obj = {f"sample{i+1}": sample for i, sample in enumerate(samples)}

temp_gene_id = ggc_df.index
gene_id = temp_gene_id.tolist()
top_gene_id_list = top_gene_list.iloc[:,0].tolist()
count_row = len(ggc_df)
num_states = 2**(len(top_gene_list))
num_actions = count_row
terminated_index = count_row # int for loop break, index that out of range

num_episodes = 2000
record = 0
total_score = 0
selected_gene_counter = []

plot_connectivity = []
plot_mean_scores = []
plot_rewards = []

# Training loop
for key, value in sample_obj.items():
    print(f"Processing {key}...")
    env = GeneEnv(ggc_df, top_gene_id_list, count_row, gene_id, value)
    agent = Agent(num_states, num_actions, num_episodes, terminated_index)
    top_gene_list_index = env.get_index_top_genes()
    
    for episodes in range(num_episodes):
        terminated = False
        log_action = [terminated_index]  # Reset log for each episode
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
                # Handle case where action is repeated
                action = terminated_index
                terminated, selection, gene_expressed = env.play_step(action, terminated_index)
                log_action.append(action)
                reward = agent.cal_reward(log_action, env, selected_gene, top_gene_list_index, terminated, gene_expressed)
                cul_reward += reward
                state_new = agent.get_state(env, terminated_index, top_gene_list_index)
                state_new_index = agent.get_state_index(state_new)
                agent.update_Qtable(state_old_index, action, reward, state_new_index, terminated)
        
        # Log cumulative reward and connectivity
        plot_rewards.append(cul_reward)
        connectivity = env.cal_max_connectivity(selection)
        if connectivity > record:
            record = connectivity
                
        plot_connectivity.append(connectivity)
        total_score += connectivity
        plot_mean_scores.append(total_score / (episodes + 1))
        
        print(f"Episode {episodes}: Cumulative Reward: {cul_reward}, Connectivity: {connectivity}, Record: {record}")

        # Track selected genes in this episode
        selected_gene_counter.append(log_action)

        # Reset the environment after each episode
        env.reset()

# After training is done, plot the results
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
# plt.show()
plt.savefig(f'training_GSE108110_results_plot_{timestamp}.png')

# Check the structure and length of each row
for i, row in enumerate(selected_gene_counter):
    print(f"Row {i} length: {len(row)}")

# You can also print the first few rows to inspect
print(selected_gene_counter[:5])

max_length = max(len(row) for row in selected_gene_counter)

# Pad rows dynamically to match the longest row (with None values for missing entries)
padded_selected_gene_counter = [
    row + [None] * (max_length - len(row)) for row in selected_gene_counter
]

# Create DataFrame with dynamic column names based on the longest row
selected_gene_df = pd.DataFrame(
    padded_selected_gene_counter, 
    columns=[f"Gene{idx+1}" for idx in range(max_length)]
)

# Save the DataFrame to a CSV file
selected_gene_df.to_csv('final_v2_cross_sample_selected_gene_count_stateGSE108110_2000.csv', index=False)

q_table = agent.model.get_Q_table()
q_table_df = pd.DataFrame(q_table)
q_table_df.to_csv('final_v2_cross_sample_q-table_2000_GSE108110.csv')

et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
