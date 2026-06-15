# -*- coding: utf-8 -*-
"""Interaktīva 2D granulārās gāzes simulācija.

Programma darbojas bez papildu projekta struktūras, ja ir instalēti numpy un
matplotlib:

    python granular_gas_simulation.py

Fizikas kodola ātrai pārbaudei bez grafiskās saskarnes:

    python granular_gas_simulation.py --self-test
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict, deque
from typing import DefaultDict, Deque, Iterator, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.widgets import Button, Slider
import numpy as np


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8")


VectorArray = np.ndarray
CellKey = Tuple[int, int]


class GranularSimulation:
    """2D granulārās gāzes fizikas kodols.

    Modelis:
    - kvadrātiska kaste ar izmēru ``box_size``;
    - kastes sienas ir absolūti elastīgas;
    - gravitācija netiek ņemta vērā;
    - visas daļiņas ir vienādas: rādiuss ``radius``, masa ``mass``;
    - daļiņu sadursmes nosaka restitūcijas koeficients ``restitution``
      diapazonā [0, 1].

    Sistēmas stāvoklis glabājas blīvos numpy masīvos:
    - ``pos``: daļiņu koordinātas, forma (N, 2);
    - ``vel``: daļiņu ātrumi, forma (N, 2).

    Sadursmju meklēšanai tiek izmantots cell list jeb telpiskais režģis. Ja
    daļiņu blīvums ir aptuveni ierobežots, katrai daļiņai ir tikai neliels
    lokālo kaimiņu skaits, tāpēc viens solis praktiski mērogojas kā O(N), nevis
    O(N^2).
    """

    def __init__(
        self,
        num_particles: int = 120,
        box_size: float = 1.0,
        radius: float = 0.014,
        mass: float = 1.0,
        restitution: float = 0.90,
        dt: float = 0.0025,
        initial_speed: float = 1.0,
        seed: Optional[int] = None,
        positions: Optional[VectorArray] = None,
        velocities: Optional[VectorArray] = None,
        cell_size: Optional[float] = None,
        collision_iterations: int = 2,
    ) -> None:
        self.n = int(num_particles)
        self.box_size = float(box_size)
        self.radius = float(radius)
        self.mass = float(mass)
        self.restitution = float(restitution)
        self.dt = float(dt)
        self.initial_speed = float(initial_speed)
        self.collision_iterations = int(collision_iterations)
        self.rng = np.random.default_rng(seed)

        self._validate_parameters()

        diameter = 2.0 * self.radius
        self.cell_size = diameter if cell_size is None else float(cell_size)
        self.cell_size = max(self.cell_size, diameter)
        self.num_cells = max(1, int(np.ceil(self.box_size / self.cell_size)))

        self.pos = self._init_positions(positions)
        self.vel = self._init_velocities(velocities)

        self.time = 0.0
        self.step_index = 0

    def step(self) -> int:
        """Izpildīt vienu integrācijas soli un atgriezt sadursmju skaitu."""

        self.pos += self.vel * self.dt
        self._resolve_wall_collisions()

        collision_count = 0
        for _ in range(self.collision_iterations):
            collision_count += self.resolve_collisions()
            self._resolve_wall_collisions()

        self.time += self.dt
        self.step_index += 1
        return collision_count

    def resolve_collisions(self) -> int:
        """Atrast un atrisināt daļiņu sadursmes, izmantojot telpisko režģi."""

        collision_count = 0
        min_distance = 2.0 * self.radius
        min_distance_sq = min_distance * min_distance

        for i, j in self._iter_candidate_pairs():
            delta_pos = self.pos[j] - self.pos[i]
            distance_sq = float(np.dot(delta_pos, delta_pos))

            if distance_sq >= min_distance_sq:
                continue

            if distance_sq > 1e-24:
                distance = float(np.sqrt(distance_sq))
                normal = delta_pos / distance
            else:
                # Ļoti reta situācija, kad centri gandrīz sakrīt. Nejaušs
                # virziens ļauj ģeometriski atdalīt daļiņas bez dalīšanas ar 0.
                distance = 0.0
                normal = self._random_unit_vector()

            self._separate_overlapping_particles(
                i=i,
                j=j,
                normal=normal,
                overlap=min_distance - distance,
            )

            relative_velocity = self.vel[j] - self.vel[i]

            # Galvenā sadursmes formula:
            #
            # n ir vienības normāle no daļiņas i centra uz daļiņas j centru.
            # Relatīvais ātrums v_rel = v_j - v_i tiek sadalīts normālajā un
            # tangenciālajā komponentē.
            #
            # Matricu formā normālo projekciju dod projektors:
            #     P_n = n n^T
            #     v_n = P_n v_rel = n (v_rel · n)
            #
            # Tangenciālā daļa (I - P_n) v_rel nemainās, jo modelī nav berzes.
            # Restitūcijas koeficients e maina tikai normālo komponenti:
            #     (v_rel' · n) = -e (v_rel · n)
            #
            # Vienādu masu gadījumā no impulsa saglabāšanās seko, ka impulsa
            # korekcija starp abām daļiņām sadalās vienādi.
            normal_velocity = float(np.dot(relative_velocity, normal))

            if normal_velocity >= 0.0:
                # Daļiņas jau attālinās normāles virzienā. Otru impulsu
                # pielietot nedrīkst, citādi enerģija mākslīgi pieaugtu.
                continue

            impulse = (
                0.5
                * (1.0 + self.restitution)
                * normal_velocity
                * normal
            )

            self.vel[i] += impulse
            self.vel[j] -= impulse
            collision_count += 1

        return collision_count

    def get_kinetic_energy(self) -> float:
        """Atgriezt sistēmas kopējo kinētisko enerģiju."""

        speed_sq = np.einsum("ij,ij->i", self.vel, self.vel)
        return float(0.5 * self.mass * np.sum(speed_sq))

    def get_total_kinetic_energy(self) -> float:
        """Saderības metode ar iepriekšējo nosaukumu."""

        return self.get_kinetic_energy()

    def shake(self, amplitude: float = 0.55) -> None:
        """Pievienot daļiņām nejaušu termisko ātruma impulsu.

        Poga "Sakratīt" modelē ārēju enerģijas pievadīšanu. No nejaušā impulsa
        tiek atņemta vidējā vērtība, lai neradītu mākslīgu masas centra dreifu.
        """

        velocity_kick = self.rng.normal(
            loc=0.0,
            scale=amplitude,
            size=self.vel.shape,
        )
        velocity_kick -= np.mean(velocity_kick, axis=0)
        self.vel += velocity_kick

    def reset_clock(self) -> None:
        """Nodzēst simulācijas pulksteni, nemainot daļiņu stāvokli."""

        self.time = 0.0
        self.step_index = 0

    def _validate_parameters(self) -> None:
        if self.n <= 0:
            raise ValueError("num_particles jābūt pozitīvam.")
        if self.box_size <= 0.0:
            raise ValueError("box_size jābūt pozitīvam.")
        if self.radius <= 0.0:
            raise ValueError("radius jābūt pozitīvam.")
        if self.mass <= 0.0:
            raise ValueError("mass jābūt pozitīvam.")
        if not 0.0 <= self.restitution <= 1.0:
            raise ValueError("restitution jābūt diapazonā [0, 1].")
        if self.dt <= 0.0:
            raise ValueError("dt jābūt pozitīvam.")
        if self.collision_iterations <= 0:
            raise ValueError("collision_iterations jābūt pozitīvam.")
        if 2.0 * self.radius >= self.box_size:
            raise ValueError("Daļiņas diametram jābūt mazākam par kasti.")

    def _init_positions(
        self,
        positions: Optional[VectorArray],
    ) -> VectorArray:
        if positions is not None:
            array = np.array(positions, dtype=float, copy=True)
            self._validate_state_array(array, "positions")
            return array

        return self._generate_non_overlapping_positions()

    def _init_velocities(
        self,
        velocities: Optional[VectorArray],
    ) -> VectorArray:
        if velocities is not None:
            array = np.array(velocities, dtype=float, copy=True)
            self._validate_state_array(array, "velocities")
            return array

        angles = self.rng.uniform(0.0, 2.0 * np.pi, size=self.n)
        velocity = self.initial_speed * np.column_stack(
            (np.cos(angles), np.sin(angles))
        )

        # Noņemam kopējo impulsu, lai gāze nelidotu kā viens mākonis.
        velocity -= np.mean(velocity, axis=0)
        rms_speed = float(np.sqrt(np.mean(np.sum(velocity * velocity, axis=1))))
        if rms_speed > 0.0:
            velocity *= self.initial_speed / rms_speed

        return velocity

    def _validate_state_array(self, array: VectorArray, name: str) -> None:
        if array.shape != (self.n, 2):
            raise ValueError(f"{name} masīvam jābūt ar formu ({self.n}, 2).")
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{name} satur NaN vai bezgalīgas vērtības.")

    def _generate_non_overlapping_positions(self) -> VectorArray:
        positions = np.empty((self.n, 2), dtype=float)
        low = self.radius
        high = self.box_size - self.radius
        min_distance_sq = (2.0 * self.radius) ** 2
        max_attempts = 40_000

        for particle_index in range(self.n):
            for _ in range(max_attempts):
                candidate = self.rng.uniform(low, high, size=2)

                if particle_index == 0:
                    positions[particle_index] = candidate
                    break

                delta = positions[:particle_index] - candidate
                distance_sq = np.einsum("ij,ij->i", delta, delta)
                if np.all(distance_sq >= min_distance_sq):
                    positions[particle_index] = candidate
                    break
            else:
                raise ValueError(
                    "Neizdevās izvietot daļiņas bez pārklāšanās. "
                    "Samaziniet N vai rādiusu R."
                )

        return positions

    def _resolve_wall_collisions(self) -> None:
        """Atrisināt absolūti elastīgas sadursmes ar kastes sienām."""

        left = self.pos[:, 0] < self.radius
        right = self.pos[:, 0] > self.box_size - self.radius
        bottom = self.pos[:, 1] < self.radius
        top = self.pos[:, 1] > self.box_size - self.radius

        self.pos[left, 0] = self.radius
        self.vel[left, 0] = np.abs(self.vel[left, 0])

        self.pos[right, 0] = self.box_size - self.radius
        self.vel[right, 0] = -np.abs(self.vel[right, 0])

        self.pos[bottom, 1] = self.radius
        self.vel[bottom, 1] = np.abs(self.vel[bottom, 1])

        self.pos[top, 1] = self.box_size - self.radius
        self.vel[top, 1] = -np.abs(self.vel[top, 1])

    def _build_cell_list(self) -> DefaultDict[CellKey, List[int]]:
        """Sadalīt daļiņas telpiskā režģa šūnās."""

        grid: DefaultDict[CellKey, List[int]] = defaultdict(list)
        cell_indices = np.floor(self.pos / self.cell_size).astype(int)
        cell_indices = np.clip(cell_indices, 0, self.num_cells - 1)

        for particle_index, (cell_x, cell_y) in enumerate(cell_indices):
            grid[(int(cell_x), int(cell_y))].append(particle_index)

        return grid

    def _iter_candidate_pairs(self) -> Iterator[Tuple[int, int]]:
        """Ģenerēt tikai pārus no vienas vai blakus esošām šūnām.

        Ja šūnas izmērs nav mazāks par daļiņas diametru, jebkurš potenciāli
        sadurošs pāris atrodas tajā pašā vai vienā no astoņām blakus šūnām.
        Tāpēc nav vajadzīgs pilns N(N-1)/2 pāru pārlasījums.
        """

        grid = self._build_cell_list()
        neighbor_offsets = (
            (-1, -1),
            (-1, 0),
            (-1, 1),
            (0, -1),
            (0, 0),
            (0, 1),
            (1, -1),
            (1, 0),
            (1, 1),
        )

        for (cell_x, cell_y), particle_indices in grid.items():
            for i in particle_indices:
                for offset_x, offset_y in neighbor_offsets:
                    neighbor_key = (cell_x + offset_x, cell_y + offset_y)
                    for j in grid.get(neighbor_key, []):
                        if j > i:
                            yield i, j

    def _separate_overlapping_particles(
        self,
        i: int,
        j: int,
        normal: VectorArray,
        overlap: float,
    ) -> None:
        """Noņemt ģeometrisku pārklāšanos pēc diskrētā laika soļa."""

        correction = 0.5 * overlap * normal
        self.pos[i] -= correction
        self.pos[j] += correction

    def _random_unit_vector(self) -> VectorArray:
        vector = self.rng.normal(size=2)
        norm = float(np.linalg.norm(vector))
        if norm == 0.0:
            return np.array([1.0, 0.0])
        return vector / norm


class GranularGasGUI:
    """Zinātniska GUI: animācija, analītika un vadības panelis."""

    @staticmethod
    def _normalize_energy_curve(
        energies: np.ndarray,
        initial_energy: float,
    ) -> np.ndarray:
        """Normalizēt kinētisko enerģiju uz sākuma līmeni katram palaidumam."""

        energies = np.asarray(energies, dtype=float)
        if initial_energy <= 0.0:
            return np.zeros_like(energies)
        return energies / float(initial_energy)

    def __init__(
        self,
        num_particles: int = 120,
        initial_restitution: float = 0.90,
        seed: int = 7,
        history_size: int = 12_000,
        interval_ms: int = 20,
        radius: float = 0.014,
    ) -> None:
        initial_restitution = float(np.clip(initial_restitution, 0.50, 1.00))

        self.simulation = GranularSimulation(
            num_particles=num_particles,
            radius=radius,
            restitution=initial_restitution,
            seed=seed,
        )
        self.history_size = history_size
        self.interval_ms = interval_ms
        self.last_collision_count = 0

        self.times: Deque[float] = deque(maxlen=history_size)
        self.energies: Deque[float] = deque(maxlen=history_size)
        self.initial_energy = self.simulation.get_kinetic_energy()
        self.saved_curve_lines: List[Line2D] = []
        self.saved_curve_count = 0

        self.figure = plt.figure(figsize=(13.2, 7.6))
        manager = self.figure.canvas.manager
        if manager is not None:
            manager.set_window_title("2D granulārā gāze")

        self._build_layout()
        self._create_artists()
        self._create_widgets(initial_restitution)
        self._connect_widgets()
        self._reset_energy_history()

        self.animation = FuncAnimation(
            self.figure,
            self._update_frame,
            interval=self.interval_ms,
            blit=False,
            cache_frame_data=False,
        )

    def show(self) -> None:
        """Parādīt lietotnes logu."""

        plt.show()

    def _build_layout(self) -> None:
        grid_spec = self.figure.add_gridspec(
            nrows=2,
            ncols=2,
            height_ratios=(5, 1),
            left=0.06,
            right=0.96,
            top=0.90,
            bottom=0.11,
            hspace=0.34,
            wspace=0.26,
        )

        self.ax_box = self.figure.add_subplot(grid_spec[0, 0])
        self.ax_energy = self.figure.add_subplot(grid_spec[0, 1])
        self.ax_controls = self.figure.add_subplot(grid_spec[1, :])
        self.ax_controls.axis("off")

        self.figure.suptitle(
            "2D granulārā gāze: telpiskais režģis, restitūcija, enerģija",
            fontsize=14,
        )

        self.ax_box.set_title("Daļiņu animācija")
        self.ax_box.set_xlim(0.0, self.simulation.box_size)
        self.ax_box.set_ylim(0.0, self.simulation.box_size)
        self.ax_box.set_aspect("equal", adjustable="box")
        self.ax_box.set_xlabel("x")
        self.ax_box.set_ylabel("y")

        self.ax_energy.set_title("Kinētiskās enerģijas evolūcija")
        self.ax_energy.set_xlabel("Laiks t (modeļa vienības)")
        self.ax_energy.set_ylabel("E_k / E_k(0) (normēta enerģija)")
        self.ax_energy.grid(True, alpha=0.30)

        self.status_text = self.ax_controls.text(
            0.02,
            0.18,
            "",
            transform=self.ax_controls.transAxes,
            fontsize=10,
            color="#303030",
        )

    def _create_artists(self) -> None:
        colors = plt.cm.viridis(np.linspace(0.12, 0.95, self.simulation.n))
        self.circles: List[Circle] = []

        for position, color in zip(self.simulation.pos, colors):
            circle = Circle(
                xy=position,
                radius=self.simulation.radius,
                facecolor=color,
                edgecolor="black",
                linewidth=0.30,
                alpha=0.92,
            )
            self.ax_box.add_patch(circle)
            self.circles.append(circle)

        self.energy_line, = self.ax_energy.plot(
            [],
            [],
            color="#c83f3f",
            linewidth=2.0,
            label="aktuālā trajektorija",
        )

    def _create_widgets(self, initial_restitution: float) -> None:
        restitution_ax = self.figure.add_axes([0.17, 0.055, 0.30, 0.035])
        speed_ax = self.figure.add_axes([0.17, 0.018, 0.30, 0.035])
        shake_ax = self.figure.add_axes([0.52, 0.044, 0.12, 0.045])
        save_ax = self.figure.add_axes([0.66, 0.044, 0.14, 0.045])
        clear_ax = self.figure.add_axes([0.82, 0.044, 0.11, 0.045])

        self.restitution_slider = Slider(
            ax=restitution_ax,
            label="restitūcija e",
            valmin=0.50,
            valmax=1.00,
            valinit=initial_restitution,
            valstep=0.01,
        )
        self.steps_per_frame_slider = Slider(
            ax=speed_ax,
            label="soļi/kadrā",
            valmin=1,
            valmax=40,
            valinit=1,
            valstep=1,
        )
        self.shake_button = Button(
            ax=shake_ax,
            label="Sakratīt",
            hovercolor="#f0e5d8",
        )
        self.save_curve_button = Button(
            ax=save_ax,
            label="Saglabāt līkni",
            hovercolor="#e7eef8",
        )
        self.clear_curves_button = Button(
            ax=clear_ax,
            label="Notīrīt",
            hovercolor="#f4dddd",
        )

    def _connect_widgets(self) -> None:
        self.restitution_slider.on_changed(self._set_restitution)
        self.shake_button.on_clicked(self._shake_particles)
        self.save_curve_button.on_clicked(self._save_current_curve)
        self.clear_curves_button.on_clicked(self._clear_saved_curves)

    def _reset_energy_history(self) -> None:
        self.times.clear()
        self.energies.clear()
        self.times.append(0.0)
        self.energies.append(self.simulation.get_kinetic_energy())
        self.initial_energy = self.energies[-1]

    def _set_restitution(self, value: float) -> None:
        self.simulation.restitution = float(value)
        self._update_status()

    def _shake_particles(self, _event: object) -> None:
        self.simulation.shake(amplitude=0.55)
        self.simulation.reset_clock()
        self._reset_energy_history()
        self._update_energy_plot()
        self._update_status()

    def _save_current_curve(self, _event: object) -> None:
        times = np.fromiter(self.times, dtype=float)
        energies = np.fromiter(self.energies, dtype=float)
        if len(times) < 2:
            return

        self.saved_curve_count += 1
        normalized_energies = self._normalize_energy_curve(
            energies,
            self.initial_energy,
        )
        label = (
            f"variants {self.saved_curve_count}: "
            f"e={self.simulation.restitution:.2f}, "
            f"E0={self.initial_energy:.2f}"
        )
        line, = self.ax_energy.plot(
            times.copy(),
            normalized_energies.copy(),
            linewidth=1.6,
            alpha=0.72,
            label=label,
        )
        self.saved_curve_lines.append(line)
        self._refresh_legend()
        self._update_energy_plot()

    def _clear_saved_curves(self, _event: object) -> None:
        for line in self.saved_curve_lines:
            line.remove()
        self.saved_curve_lines.clear()
        self.saved_curve_count = 0
        legend = self.ax_energy.get_legend()
        if legend is not None:
            legend.remove()
        self._update_energy_plot()

    def _update_frame(self, _frame: int) -> List[object]:
        self.simulation.restitution = float(self.restitution_slider.val)
        steps_per_frame = int(self.steps_per_frame_slider.val)
        collision_count = 0

        for _ in range(steps_per_frame):
            collision_count += self.simulation.step()
            self.times.append(self.simulation.time)
            self.energies.append(self.simulation.get_kinetic_energy())

        self.last_collision_count = collision_count

        for circle, position in zip(self.circles, self.simulation.pos):
            circle.center = position

        self._update_energy_plot()
        self._update_status()
        return [*self.circles, self.energy_line, self.status_text]

    def _update_energy_plot(self) -> None:
        times = np.fromiter(self.times, dtype=float)
        energies = np.fromiter(self.energies, dtype=float)
        normalized_energies = self._normalize_energy_curve(
            energies,
            self.initial_energy,
        )
        self.energy_line.set_data(times, normalized_energies)

        time_min = 0.0
        time_max = max(float(times[-1]), 0.1)
        energy_max = float(np.max(normalized_energies)) if normalized_energies.size else 1.0

        for line in self.saved_curve_lines:
            saved_x = np.asarray(line.get_xdata(), dtype=float)
            saved_y = np.asarray(line.get_ydata(), dtype=float)
            if saved_x.size:
                time_max = max(time_max, float(np.max(saved_x)))
            if saved_y.size:
                energy_max = max(energy_max, float(np.max(saved_y)))

        self.ax_energy.set_xlim(time_min, time_max)
        self.ax_energy.set_ylim(0.0, max(energy_max * 1.08, 1.0))

    def _refresh_legend(self) -> None:
        if not self.saved_curve_lines:
            return

        self.ax_energy.legend(
            loc="upper right",
            fontsize=8,
            framealpha=0.88,
        )

    def _update_status(self) -> None:
        energy = self.energies[-1]
        energy_ratio = energy / self.initial_energy if self.initial_energy else 0.0
        self.status_text.set_text(
            f"N = {self.simulation.n}; "
            f"R = {self.simulation.radius:.3f}; "
            f"e = {self.simulation.restitution:.2f}; "
            f"t = {self.simulation.time:.3f}; "
            f"E/E0 = {energy_ratio:.4f}; "
            f"soļi/kadrā = {int(self.steps_per_frame_slider.val)}; "
            f"sadursmes/kadrā = {self.last_collision_count}; "
            f"saglabātas līknes = {len(self.saved_curve_lines)}"
        )


def run_self_test() -> None:
    """Minimālas fizikas pārbaudes bez pytest."""

    normalized = GranularGasGUI._normalize_energy_curve(
        np.array([2.0, 4.0, 6.0], dtype=float),
        initial_energy=2.0,
    )
    np.testing.assert_allclose(normalized, [1.0, 2.0, 3.0])

    elastic = GranularSimulation(
        num_particles=2,
        box_size=1.0,
        radius=0.05,
        restitution=1.0,
        positions=np.array([[0.45, 0.50], [0.54, 0.50]]),
        velocities=np.array([[1.0, 0.0], [0.0, 0.0]]),
    )

    energy_before = elastic.get_kinetic_energy()
    elastic.resolve_collisions()
    energy_after = elastic.get_kinetic_energy()

    np.testing.assert_allclose(elastic.vel[0], [0.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(elastic.vel[1], [1.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(energy_before, energy_after, rtol=1e-12)

    inelastic = GranularSimulation(
        num_particles=2,
        box_size=1.0,
        radius=0.05,
        restitution=0.8,
        positions=np.array([[0.45, 0.50], [0.54, 0.50]]),
        velocities=np.array([[1.0, 0.0], [0.0, 0.0]]),
    )

    energy_before = inelastic.get_kinetic_energy()
    inelastic.resolve_collisions()
    energy_after = inelastic.get_kinetic_energy()

    if not energy_after < energy_before:
        raise AssertionError("Neelastīgai sadursmei jāsamazina enerģija.")

    clock_test = GranularSimulation(num_particles=20, restitution=0.9, seed=11)
    for _ in range(10):
        clock_test.step()
    clock_test.shake()
    clock_test.reset_clock()
    np.testing.assert_allclose(clock_test.time, 0.0)

    for restitution in (1.0, 0.9):
        simulation = GranularSimulation(
            num_particles=50,
            radius=0.013,
            restitution=restitution,
            seed=3,
        )
        energy_start = simulation.get_kinetic_energy()
        for _ in range(250):
            simulation.step()
        energy_end = simulation.get_kinetic_energy()
        print(
            f"restitūcija={restitution:.1f}: "
            f"E0={energy_start:.8f}, "
            f"E250={energy_end:.8f}, "
            f"attiecība={energy_end / energy_start:.8f}"
        )

    print("pašpārbaude pabeigta")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Interaktīva 2D granulārās gāzes simulācija."
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Pārbaudīt fizikas kodolu bez GUI atvēršanas.",
    )
    parser.add_argument(
        "--particles",
        type=int,
        default=120,
        help="Sākotnējais daļiņu skaits GUI logā.",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=0.014,
        help="Sākotnējais daļiņas rādiuss GUI logā.",
    )
    parser.add_argument(
        "--restitution",
        type=float,
        default=0.90,
        help="Sākotnējais restitūcijas koeficients GUI logā.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.self_test:
        run_self_test()
        return

    app = GranularGasGUI(
        num_particles=args.particles,
        radius=args.radius,
        initial_restitution=args.restitution,
        seed=7,
    )
    app.show()


if __name__ == "__main__":
    main()
