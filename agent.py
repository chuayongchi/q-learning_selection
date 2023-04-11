import random
import numpy as np
from model import Q_learning

NUMBER_GENE = 20

class Agent():

    def __init__(self, num_states, num_actions, num_episodes, terminated_index):
        self.epsilon = 0.2 # randomness
        self.model = Q_learning(num_states, num_actions, num_episodes)
        self.num_actions = num_actions
        self.selected_gene = terminated_index + 1
        self.previous_selected_gene = 0

    

    def get_action(self, state):
        ## AGENT NEED TO DECIDE WHICH SLOT HAVE TO WRITE 1
        # initialize the action array with the size of the number of gene
        # set the action randomlly select the gene
        # random moves: exploration / 
        # best moves (based on Q-value): exploitation
            
        if random.uniform(0, 1) < self.epsilon:
            action = np.random.randint(0,self.num_actions)
            print(action)
        else:
            possibleQs = self.model.Q[state, :]
            action = np.argmax(possibleQs)
        return action

    def get_state(self, env, terminated_index, top_gene_list_index):
        previous_selected_gene = self.previous_selected_gene
        if (self.selected_gene == (terminated_index + 1)):
            state = [(False),(False),(False),(False),(False),(False),(False),(False),(False),(False),(False)]
        
        elif(self.selected_gene == terminated_index):
            state = [
                (env.ggc[top_gene_list_index[0]][previous_selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[1]][previous_selected_gene] > 0.1),
                
                (env.ggc[top_gene_list_index[2]][previous_selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[3]][previous_selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[4]][previous_selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[5]][previous_selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[6]][previous_selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[7]][previous_selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[8]][previous_selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[9]][previous_selected_gene] > 0.1)
            ]

        else:
            state = [
                (env.ggc[top_gene_list_index[0]][self.selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[1]][self.selected_gene] > 0.1),
                
                (env.ggc[top_gene_list_index[2]][self.selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[3]][self.selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[4]][self.selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[5]][self.selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[6]][self.selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[7]][self.selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[8]][self.selected_gene] > 0.1),

                (env.ggc[top_gene_list_index[9]][self.selected_gene] > 0.1)
            ]
            self.previous_selected_gene = self.selected_gene

        return np.array(state, dtype=int)

    def get_state_index(self, state):
        state_index = 0
        for i in range(10):
            state_index += 2**i*state[i]
        return state_index
    
    def cal_reward(self,action_log, env, selected_gene, top_gene_list_index, terminated, gene_expressed):
        reward =0
        log = action_log

        m = len(top_gene_list_index)
        n = len(log)-1
        print("gene_expressed:", gene_expressed)
        print("top gene list:", top_gene_list_index)

        if (not terminated and (gene_expressed == 1)):
            print("n:", n)
            for j in range(m):
                if (env.ggc[top_gene_list_index[j]][selected_gene] < 0.01):
                    reward = reward + (-1.5)*(env.ggc[top_gene_list_index[j]][selected_gene])
                else:
                    print("selected gene:", selected_gene)
                    print("top gene:", top_gene_list_index[j] )
                    reward = reward + (env.ggc[top_gene_list_index[j]][selected_gene])
                print("inside agent, the cal rewards, the reward is ", reward)
            
        return reward
    
    def get_selected_gene(self, action):
        self.selected_gene = action
        return self.selected_gene

    def update_Qtable(self, state_old, action, reward, state_new, terminated):
        self.model.update_Qvalue(state_old, action, reward, state_new, terminated)

    def update_Qs_table(self, num_episodes):
        self.model.update_Qs_value(num_episodes)
    
    