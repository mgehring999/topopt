import gymnasium as gym
import torch
import torch.nn as nn
import numpy as np
import time

from topopt.rl import DQNNetwork

# Load network and weights
obs_dim = 4
n_actions = 2
q_net = DQNNetwork(obs_dim, n_actions)
q_net.load_state_dict(torch.load("dqn_cartpole.pth"))
q_net.eval()

env = gym.make("CartPole-v1", render_mode="human")

for ep in range(5):
    obs = env.reset()[0]
    done = False
    total_reward = 0

    while not done:
        obs_tensor = torch.tensor(obs, dtype=torch.float32).unsqueeze(0)
        with torch.no_grad():
            action = q_net(obs_tensor).argmax().item()
        obs, reward, done, _, _ = env.step(action)
        total_reward += reward
        time.sleep(0.02)

    print(f"Test episode {ep+1}: Total reward = {total_reward}")

env.close()
