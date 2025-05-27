# modius-strip
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class MobiusStrip:
    def __init__(self, R=1.0, w=0.2, n=100):
        self.R = R
        self.w = w
        self.n = n
        self.u, self.v = np.meshgrid(
            np.linspace(0, 2 * np.pi, n),
            np.linspace(-w / 2, w / 2, n)
        )
        self.x, self.y, self.z = self._compute_coordinates()

    def _compute_coordinates(self):
        u = self.u
        v = self.v
        R = self.R
        x = (R + v * np.cos(u / 2)) * np.cos(u)
        y = (R + v * np.cos(u / 2)) * np.sin(u)
        z = v * np.sin(u / 2)
        return x, y, z

    def compute_surface_area(self):
        dx_u = np.gradient(self.x, axis=0)
        dx_v = np.gradient(self.x, axis=1)
        dy_u = np.gradient(self.y, axis=0)
        dy_v = np.gradient(self.y, axis=1)
        dz_u = np.gradient(self.z, axis=0)
        dz_v = np.gradient(self.z, axis=1)

      
        cross_x = dy_u * dz_v - dz_u * dy_v
        cross_y = dz_u * dx_v - dx_u * dz_v
        cross_z = dx_u * dy_v - dy_u * dx_v

        dA = np.sqrt(cross_x**2 + cross_y**2 + cross_z**2)
        area = np.sum(dA) * (2 * np.pi / self.n) * (self.w / self.n)
        return area

    def compute_edge_length(self):
       
        u = np.linspace(0, 2 * np.pi, self.n)
        edge_lengths = []
        for v_edge in [-self.w/2, self.w/2]:
            x = (self.R + v_edge * np.cos(u / 2)) * np.cos(u)
            y = (self.R + v_edge * np.cos(u / 2)) * np.sin(u)
            z = v_edge * np.sin(u / 2)
            dx = np.gradient(x)
            dy = np.gradient(y)
            dz = np.gradient(z)
            ds = np.sqrt(dx**2 + dy**2 + dz**2)
            edge_lengths.append(np.sum(ds))
        return np.mean(edge_lengths)

    def plot(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(self.x, self.y, self.z, rstride=1, cstride=1, color='lightblue', edgecolor='k', alpha=0.8)
        ax.set_title('Möbius Strip')
        plt.show()


strip = MobiusStrip(R=1.0, w=0.4, n=200)
print("Surface Area ≈", strip.compute_surface_area())
print("Edge Length ≈", strip.compute_edge_length())
strip.plot()
