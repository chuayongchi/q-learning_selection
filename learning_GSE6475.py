import pandas as pd
import numpy as np
from env import GeneEnv
from agent import Agent
import matplotlib.pyplot as plt

# Load the data
df = pd.read_csv('./data/data_degenes_with_acne_sample.csv', index_col=0)
ggc_df = pd.read_csv('./data/WGCNA_Matrix_GSE6475_Final.csv', index_col=0)
top_gene_list = pd.read_csv('./data/gene_sign_list_GSE6475.csv')

# Prepare the sample data
sample1 = np.array(df.iloc[0].values.tolist())
sample2 = np.array(df.iloc[1].values.tolist())
sample3 = np.array(df.iloc[2].values.tolist())
sample4 = np.array(df.iloc[3].values.tolist())
sample5 = np.array(df.iloc[4].values.tolist())
sample6 = np.array(df.iloc[5].values.tolist())

sample_obj = {
    "sample1": sample1,
    "sample2": sample2,
    "sample3": sample3,
    "sample4": sample4,
    "sample5": sample5,
    "sample6": sample6,
}

# Get gene information
temp_gene_id = ggc_df.index
gene_id = temp_gene_id.tolist()

temp_top_gene_id = top_gene_list.iloc[:, 0]
top_gene_id_list = temp_top_gene_id.tolist()

# Environment setup
count_row = len(ggc_df)
size = len(top_gene_list)

num_states = 2 ** size
num_actions = count_row
actions_space = num_actions + 1
terminated_index = count_row  # Index that signifies termination

# Hyperparameters
num_episodes = 200
epsilon_start = 1.0
epsilon_decay = 0.995
epsilon_min = 0.01

# Initialize Q-table and logs
q_table = np.zeros((num_states, actions_space))
plot_connectivity = []
plot_mean_scores = []
plot_rewards = []

record = 0
total_score = 0
num_iteration = 0

# Main RL loop
for key, value in sample_obj.items():
    print(f"Processing {key}")
    env = GeneEnv(ggc_df, top_gene_id_list, count_row, gene_id, value)
    agent = Agent(num_states, num_actions, num_episodes, terminated_index, epsilon_start, epsilon_decay, epsilon_min)
    # agent = Agent(num_states, num_actions, num_episodes, terminated_index)

    top_gene_list_index = env.get_index_top_genes()

    for episode in range(num_episodes):
        num_iteration += 1
        cul_reward = 0
        terminated = False
        log_action = [terminated_index]  # Reset action log at the start of each episode
        selected_gene_list = []  # Track genes selected in this episode

        while not terminated:
            state_old = agent.get_state(env, terminated_index, top_gene_list_index)
            state_old_index = agent.get_state_index(state_old)

            # Use epsilon-greedy strategy to get the action
            action = agent.get_action(state_old_index)

            if action not in log_action:
                terminated, selection, gene_expressed = env.play_step(action, terminated_index)
                log_action.append(action)

                selected_gene = agent.get_selected_gene(action)
                selected_gene_list.append(selected_gene)

                if gene_expressed == 0:
                    # Gene not expressed; select a new action
                    print(f"Gene {selected_gene} is not expressed. Selecting a new action.")
                    continue  # Go back to the start of the loop

                reward = agent.cal_reward(log_action, env, selected_gene, top_gene_list_index, terminated, gene_expressed)
                cul_reward += reward

                state_new = agent.get_state(env, terminated_index, top_gene_list_index)
                state_new_index = agent.get_state_index(state_new)

                # Update Q-table
                agent.update_Qtable(state_old_index, action, reward, state_new_index, terminated)
            else:
                # Handle repeated actions
                action = terminated_index
                terminated, selection, gene_expressed = env.play_step(action, terminated_index)
                log_action.append(action)

                selected_gene = agent.get_selected_gene(action)
                selected_gene_list.append(selected_gene)

                reward = agent.cal_reward(log_action, env, selected_gene, top_gene_list_index, terminated, gene_expressed)
                cul_reward += reward

                state_new = agent.get_state(env, terminated_index, top_gene_list_index)
                state_new_index = agent.get_state_index(state_new)

                agent.update_Qtable(state_old_index, action, reward, state_new_index, terminated)
                plot_mean_scores.append(reward)

        # Store rewards and connectivity data
        plot_rewards.append(cul_reward)
        connectivity = env.cal_max_connectivity(selection)
        if connectivity > record:
            record = connectivity
        plot_connectivity.append(cul_reward)

        total_score += connectivity
        mean_score = total_score / (episode + 1)

        # Plot the progress
        plt.plot(plot_connectivity, label='Cumulative Connectivity')
        plt.plot(plot_mean_scores, label='Mean Scores')
        plt.legend()
        plt.show()

        # Decay epsilon after each episode
        agent.decay_epsilon()

        # Reset environment for next episode
        env.reset()

# Export the results
selected_gene_lists = [gene_id[i] for i in selected_gene_list if i != terminated_index]
gene_markers_count = {gene: selected_gene_lists.count(gene) for gene in set(selected_gene_lists)}

# Save final selected genes and their counts
selected_gene_count_df = pd.DataFrame(list(gene_markers_count.items()), columns=['Gene', 'Count'])
selected_gene_count_df.to_csv('final_State_cross_sample_selected_gene_count_state_1200genes.csv', index=False)

# Save Q-table
q_table_df = pd.DataFrame(agent.Q)
q_table_df.to_csv('final_State_cross_sample_Reward_q-table_1200genes.csv', index=False)

print("Training completed and files saved.")
