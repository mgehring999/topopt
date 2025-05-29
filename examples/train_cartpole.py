from topopt.rl import DQNAgent
import gymnasium as gym
import random, torch
import numpy as np

np.random.seed(42)
random.seed(42)
torch.manual_seed(42)

episodes = 2000
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

    print(f"Ep {episode+1}: Reward = {total_reward}, Steps = {env_steps}, Epsilon = {agent.epsilon:.3f}")

torch.save(agent.model.state_dict(), "dqn_cartpole.pth")
