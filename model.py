import numpy as np

class Q_learning():

    def __init__(self,num_states, num_actions, num_episodes):
        actions_space = num_actions + 1
        self.Q = np.zeros((num_states, actions_space)) ## while not terminated, this Q need to be updated
        self.gamma = 0.8
        self.learning_rate = 0.9

    def update_Qs_value(self, episodes):
        self.Qs[episodes, :, :] = np.copy(self.Q)

    
    def update_Qvalue(self, state, action, reward, next_state, terminated):
        self.Q[state, action] = self.Q[state, action] + self.learning_rate * (reward)
        print("self.Q[state, action] :", self.Q[state, action])
        if not(terminated):
            self.Q[state, action] = self.Q[state, action] + self.learning_rate * (reward + self.gamma * np.max(self.Q[next_state, :]) - self.Q[state, action])
            print("if not terminated; self.Q[state, action] :", self.Q[state, action])

    def get_Qs_table(self):
        Qs = self.Qs
        return Qs

    def get_Q_table(self):
        Q = self.Q
        return Q

    def get_best_Q(self, Q):
        best_Q = np.copy(Q)
        return best_Q
