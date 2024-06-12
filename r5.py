import math
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


AMU = 1.66e-27
CHARGE_unit = 1.6e-19
Q_alpha = 2.0 * CHARGE_unit
Q_target = 79.0 * CHARGE_unit
M_alpha = 4.0 * AMU
M_target = 197.0 * AMU
CK = 8.99e9  # 쿨롱 상수 (N·m²/C²)
dT = 1.0e-23  # 시간 간격
MAX_time = 4.0e-20  
V0 = 1.0e7  # 초기 속도

# 금 원자핵 반지름
BASE_RT = pow(197, 0.33333) * 1.0e-15
BB = BASE_RT / 2.0

def distance(a, b, c, d):
    return math.sqrt((a - c) * (a - c) + (b - d) * (b - d))

def acc(e, f, g, h):
    return e * (f - g) / (h * h * h)

def simulate_and_plot(rt_scale):
    
    trajectories = []

    cc = CK * Q_target * Q_alpha / M_alpha
    rt = BASE_RT * rt_scale

    # y0 값
    y0_values = [BB * i for i in range(1, 21, 2)]

    for y0 in y0_values:
        xa = -2e-13
        ya = y0
        vax = V0
        vay = 0.0
        xt = yt = vtx = vty = 0.0

        rr = distance(xa, ya, xt, yt)
        ax = acc(cc, xa, xt, rr)
        ay = acc(cc, ya, yt, rr)
        atx = -ax * M_alpha / M_target
        aty = -ay * M_alpha / M_target

        vax = vax + ax * dT 
        vay = vay + ay * dT 
        vtx = vtx + atx * dT 
        vty = vty + aty * dT

        trajectory = []

        time = 0.0
        output_step = 0

        while time < MAX_time:
            xa = xa + dT * vax
            ya = ya + dT * vay

            if (output_step % 30) == 0:
                trajectory.append((xa, ya))

            xt = xt + dT * vtx
            yt = yt + dT * vty
            rr = distance(xa, ya, xt, yt)

            if rr > rt:
                ax = acc(cc, xa, xt, rr)
                ay = acc(cc, ya, yt, rr)
            else:
                ax = acc(cc, xa, xt, rt)
                ay = acc(cc, ya, yt, rt)

            atx = -ax * M_alpha / M_target
            aty = -ay * M_alpha / M_target
            vax = vax + dT * ax
            vay = vay + dT * ay
            vtx = vtx + dT * atx
            vty = vty + dT * aty

            time = time + dT
            output_step += 1
        
        trajectories.append(trajectory)

    print("Simulation complete. Plotting the results.")

    #그래프 그리기
    fig, ax = plt.subplots()
    for trajectory in trajectories:
        x_vals = [point[0] for point in trajectory]
        y_vals = [point[1] for point in trajectory]
        ax.plot(x_vals, y_vals)

    # 금 원자핵의 위치를 점으로 표시, 크기 rt 값에 비례
    nucleus_size = rt * 1e16  
    ax.scatter([0], [0], color='red', s=nucleus_size, label='Gold Nucleus')

    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_title('Alpha Particle Trajectories')
    ax.legend()
    ax.grid(True)

    plt.show()


fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.25)
simulate_and_plot(1)


axcolor = 'lightgoldenrodyellow'
ax_slider = plt.axes([0.1, 0.1, 0.8, 0.03], facecolor=axcolor)
slider = Slider(ax_slider, 'rt Scale', 0.1, 30.0, valinit=1)


def update(val):
    ax.clear()
    simulate_and_plot(slider.val)

slider.on_changed(update)

plt.show()
