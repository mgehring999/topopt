import pytest
from topopt.rl import DQNAgent
import gymnasium as gym
import random, torch
import numpy as np

def test_cp_env_steps():
    np.random.seed(42)
    random.seed(42)
    torch.manual_seed(42)

    episodes = 100
    batch_size = 64

    env = gym.make("CartPole-v1")
    env.action_space.seed(42)
    obs_dim = env.observation_space.shape[0]
    n_actions = env.action_space.n

    agent = DQNAgent(obs_dim,n_actions)

    for episode in range(episodes):
        state,_ = env.reset(seed=42)
        done = False
        env_steps = 0
        total_reward = 0

        while not done:
            action = agent.select_action(state)
            next_state,reward,terminated,truncated,_ = env.step(action)
            done = terminated or truncated
            agent.remember(state,action,reward,next_state,done)

            state = next_state
            env_steps += 1
            total_reward += reward

            agent.replay(batch_size)
    
    assert env_steps > 100 # must be over 100 at episode 100
    assert np.isclose(agent.epsilon,0.099943)