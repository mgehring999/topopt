import sys
import json
import numpy as np
from PyQt5 import QtWidgets, QtGui, QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class RewardPlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure(figsize=(5, 2))
        self.ax = fig.add_subplot(111)
        super().__init__(fig)
        self.setParent(parent)

    def plot_rewards(self, rewards, current_step):
        self.ax.clear()
        self.ax.plot(rewards, label='Reward')
        self.ax.axvline(x=current_step, color='r', linestyle='--', label='Current Step')
        self.ax.legend()
        self.draw()

class RLViewer(QtWidgets.QWidget):
    def __init__(self, data_dict):
        super().__init__()
        self.data = data_dict
        self.episodes = list(self.data.keys())
        self.current_episode_idx = 0
        self.current_step = 0

        self.init_ui()
        self.load_episode()

    def init_ui(self):
        layout = QtWidgets.QVBoxLayout(self)

        # Episode label
        self.episode_label = QtWidgets.QLabel(self)
        self.episode_label.setAlignment(QtCore.Qt.AlignCenter)
        font = self.episode_label.font()
        font.setPointSize(14)
        self.episode_label.setFont(font)
        layout.addWidget(self.episode_label)

        # Image display
        self.image_label = QtWidgets.QLabel(self)
        self.image_label.setFixedSize(400, 400)
        self.image_label.setAlignment(QtCore.Qt.AlignCenter)
        image_container = QtWidgets.QWidget()
        image_container_layout = QtWidgets.QVBoxLayout(image_container)
        image_container_layout.addWidget(self.image_label, alignment=QtCore.Qt.AlignCenter)
        layout.addWidget(image_container)

        # Reward plot
        self.reward_canvas = RewardPlotCanvas(self)
        layout.addWidget(self.reward_canvas)

        # Navigation buttons
        nav_layout = QtWidgets.QHBoxLayout()
        self.prev_btn = QtWidgets.QPushButton("Previous Episode")
        self.next_btn = QtWidgets.QPushButton("Next Episode")
        nav_layout.addWidget(self.prev_btn)
        nav_layout.addWidget(self.next_btn)
        layout.addLayout(nav_layout)

        # Animation control buttons
        anim_layout = QtWidgets.QHBoxLayout()
        self.start_btn = QtWidgets.QPushButton("Start Animation")
        self.stop_btn = QtWidgets.QPushButton("Stop Animation")
        anim_layout.addWidget(self.start_btn)
        anim_layout.addWidget(self.stop_btn)
        layout.addLayout(anim_layout)

        # Timer to step through episode
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.next_step)
        self.timer.start(500)  # update every 500 ms

        self.prev_btn.clicked.connect(self.prev_episode)
        self.next_btn.clicked.connect(self.next_episode)
        self.start_btn.clicked.connect(self.start_animation)
        self.stop_btn.clicked.connect(self.stop_animation)

    def load_episode(self):
        episode_key = self.episodes[self.current_episode_idx]
        self.episode_label.setText(f"Episode: {episode_key}")
        episode_data = self.data[episode_key]
        self.steps = episode_data
        self.rewards = [step[2] for step in self.steps]
        self.current_step = 0
        self.update_display()

    def update_display(self):
        step = self.steps[self.current_step]
        element_states = step[0][0]  # Take element_states from observation structure

        # Convert list of 0/1 to numpy array
        size = int(len(element_states) ** 0.5)
        img_array = np.array(element_states, dtype=np.uint8).reshape((size, size)) * 255
        img_array = np.flipud(img_array)  # Flip vertically to correct upside-down issue

        height, width = img_array.shape
        bytes_per_line = width
        q_img = QtGui.QImage(img_array.tobytes(), width, height, bytes_per_line, QtGui.QImage.Format_Grayscale8)
        pixmap = QtGui.QPixmap.fromImage(q_img)

        self.image_label.setPixmap(pixmap.scaled(self.image_label.size(), QtCore.Qt.KeepAspectRatio))

        self.reward_canvas.plot_rewards(self.rewards, self.current_step)

    def next_step(self):
        self.current_step += 1
        if self.current_step >= len(self.steps):
            self.current_step = 0
        self.update_display()

    def prev_episode(self):
        self.current_episode_idx = (self.current_episode_idx - 1) % len(self.episodes)
        self.load_episode()

    def next_episode(self):
        self.current_episode_idx = (self.current_episode_idx + 1) % len(self.episodes)
        self.load_episode()

    def start_animation(self):
        self.timer.start(500)

    def stop_animation(self):
        self.timer.stop()

if __name__ == '__main__':
    import os
    app = QtWidgets.QApplication(sys.argv)
    with open("trajectories.json", "r") as f:
        data = json.load(f)
    viewer = RLViewer(data)
    viewer.show()
    sys.exit(app.exec_())
