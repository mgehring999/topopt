import torch
import torch.nn as nn
import torch.optim as optim
import random, sys
import numpy as np
from collections import deque

class DQNNetwork(nn.Module):
    def __init__(self, input_size, output_size):
        super(DQNNetwork, self).__init__()
        self.fc1 = nn.Linear(input_size, 128)
        self.fc3 = nn.Linear(128, output_size)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        return self.fc3(x)

# Experience replay
class ReplayBuffer:
    def __init__(self, capacity,device="cpu"):
        self.buffer = deque(maxlen=capacity)
        self.device = device

    def push(self, transition):
        self.buffer.append(transition)

    def sample(self, batch_size):
        batch = random.sample(self.buffer, batch_size)
        s, a, r, s_next, done = zip(*batch)

        return (
            torch.tensor(s, dtype=torch.float32, device=self.device),
            torch.tensor(a, dtype=torch.int64, device=self.device).unsqueeze(1),
            torch.tensor(r, dtype=torch.float32, device=self.device).unsqueeze(1),
            torch.tensor(s_next, dtype=torch.float32, device=self.device),
            torch.tensor(done, dtype=torch.float32, device=self.device).unsqueeze(1)
        )

    def __len__(self):
        return len(self.buffer)

class DQNAgent:
    def __init__(self, state_size, action_size):
        self.state_size = state_size
        self.action_size = action_size
        self.gamma = 0.99  # discount rate
        self.epsilon = 1.  # exploration rate
        self.epsilon_decay = 0.999
        self.epsilon_min = 0.1
        self.model = DQNNetwork(state_size, action_size)
        self.target_model = DQNNetwork(state_size, action_size)
        self.optimizer = optim.Adam(self.model.parameters(), lr=0.001)

        self.target_model.load_state_dict(self.model.state_dict())
        self.target_update_counter = 0
        self.target_update_interval = 10

        max_len = 10000
        self.memory = ReplayBuffer(max_len)

    def select_action(self, state):
        if np.random.rand() <= self.epsilon:
            return np.random.choice(self.action_size)
        
        with torch.no_grad():
            state = torch.FloatTensor(state).unsqueeze(0)
            q_values = self.model(state)
            return int(torch.argmax(q_values))

    def remember(self, state, action, reward, next_state, done):
        self.memory.push((state, action, reward, next_state, done))

        # count episodes to update target network
        # without passing the episode number to the agent
        if done: self.target_update_counter += 1

    def replay(self, batch_size):
        if len(self.memory) < batch_size:
            return

        batch = self.memory.sample(batch_size)
        states, actions, rewards, next_states, dones = batch

        current_q_values = self.model(states).gather(1, actions)
        next_q_values = self.target_model(next_states).max(1,keepdim=True)[0].detach()

        target_q_values = rewards + (1 - dones) * self.gamma * next_q_values

        loss = nn.MSELoss()(current_q_values, target_q_values)
        self.loss_val = loss.item()

        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()

        if self.epsilon >= self.epsilon_min:
            self.epsilon *= self.epsilon_decay

        if self.target_update_counter==self.target_update_interval:
            self.target_model.load_state_dict(self.model.state_dict())
            self.target_update_counter = 0