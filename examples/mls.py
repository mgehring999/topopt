#!/usr/bin/env python3
import random
import torch
import json

import numpy as np
from topopt.physical import Material
from topopt.mesh import Mesh, Displacement, Force
from fem.fem import FEModel,StructuralElement
import sys, logging, glob,os
from timeit import default_timer as timer

from utils.post import plot_dof
from topopt.env import TopoEnv
from topopt.rl import DQNAgent

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from torch.utils.tensorboard import SummaryWriter

np.random.seed(42)
random.seed(42)
torch.manual_seed(42)

logger = logging.getLogger('topopt')
tb_writer= SummaryWriter()

# Mesh and Model init
ndiv = 6
mesh = Mesh()
mesh.rect_mesh(ndiv)

mat = Material()
mat.set_structural_params(2.1e5,0.3)

support = Displacement(mesh)
support.add_by_plane([1,0],-1,0)

load = Force(mesh)
load.add_by_point((1,0),(0,-100))

fem = FEModel(mesh,mat,StructuralElement)

# RL Training
episodes = 5000
batch_size = 100

env = TopoEnv(fem,support,load)
state_size = 4*ndiv**2
action_size = ndiv**2 
agent = DQNAgent(state_size,action_size)

fig,axs=plt.subplots(2,1)
colors = ["white", "grey","grey","blue"]
nodes = [0.0, 0.4, 0.6,1.0]
cmap = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))

trajectories = {}

logger.info("starts training loop")
for episode in range(episodes):
    state = env.reset()
    total_reward = 0
    reward_series=[]
    start_time = timer()
    done = False 

    while not done:
        action = agent.select_action(state)
        next_state, reward, done = env.step(action)
        agent.remember(state, action, reward, next_state, done)

        state = next_state
        total_reward += reward
        reward_series.append(total_reward)

        agent.replay(batch_size)

        if "plot" in sys.argv: 
            np.save("out.npy",env.elem_state)

        if "save" in sys.argv:
            axs[0].imshow(env.elem_state.reshape((ndiv,ndiv)),cmap=cmap,origin="lower")
            axs[1].plot(reward_series)
            plt.savefig("output/elem_state_{}.png".format(episode+1))
            for a in axs: a.cla()

        if done:
            end_time = timer()
            print(f"Episode: {episode + 1}, Total Reward: {total_reward}, Time: {end_time-start_time}, Iterations: {env.count}, Exploration rate: {agent.epsilon}")

    # write trajectories
    trajectories[episode+1] = env.history
    with open("trajectories.json","w+") as f:
        json.dump(trajectories,f)

    tb_writer.add_scalar("Total Reward",total_reward,episode)
    tb_writer.add_scalar("Number of Iterations",env.count,episode)
    tb_writer.add_scalar("Strain Reward",env.strain_reward,episode)
    tb_writer.add_scalar("Volume Reaward",env.vol_reward,episode)
    tb_writer.add_scalar("Deformation Minimum",env.umin,episode)

tb_writer.close()
