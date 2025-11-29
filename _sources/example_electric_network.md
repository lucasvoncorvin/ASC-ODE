# Electric Network

We want to model an electric network with an ODE. We are given a Voltage source $U_0 = cos(100 \pi t)$. We are given the variables $R$ for the resistance and $C$ for the capacity. 
The system can be written as :

$$
\begin{aligned}
U_c(t) + RC\frac{dU_c}{dt} &= U_0(t), \\
\frac{dU_c}{dt} &= -\frac{U_c(t)}{RC} + \frac{U_0(t)}{RC}.
\end{aligned}
$$


So we can see our system is decomposed into a decaying exponential solution and oscillatory solution. 
The decay of the exponential solution is hereby determined by the term $\frac{-1}{RC}$. If $RC$ becomes small the transient behaviour vanishes very quickly.
We can solve the system and compute exact solution by splitting into homogeneous and particular solution and define $\tau = RC$ and $\omega = 100 pi$.

$$
U_C(t)=-\frac{1}{1+(\tau\omega)^2} e^{-t/\tau}+\frac{1}{1+(\tau\omega)^2}\cos(\omega t)+\frac{\tau\omega}{1+(\tau\omega)^2}\sin(\omega t).
$$

In the first test case, we set $R=C=1$. Therefore we expect visible decay of the amplitude of the oscillations and we expect similar results for all methods, where explicit might overshoot and implicit might undershoot amplitudes.
Since our stationary frequency is given by $f = 100 \pi / 2 pi = 50 Hz$, our period $T = 0.02s$, our time step should be small enough to detect oscillatory behavior. We run the simulation from $t_0 = 0$ to $t_{end} = 5s$.
So we expect for $N \geq t_{total}/\Delta t = 5/0.02 = 250$, that we see some oscillatory behavior.
And in general we expect to see mostly a sinusoidal function since the prefactor of the exponential and the cosinus are two orders lower then the prefactor of the sine term.

5000 time steps
```{image} ../demos/plots_rc/rc_methods_time_5000_R_C_1.png
:width: 45%
```
```{image} ../demos/plots_rc/rc_methods_phase_5000_R_C_1.png
:width: 45%
```

1000 time steps
```{image} ../demos/plots_rc/rc_methods_time_5_1000_R_C_1.png
:width: 45%
```
```{image} ../demos/plots_rc/rc_methods_phase_5_1000_R_C_1.png
:width: 45%
```

200 time steps
```{image} ../demos/plots_rc/rc_methods_time_5_200_R_C_1.png
:width: 45%
```
```{image} ../demos/plots_rc/rc_methods_phase_5_200_R_C_1.png
:width: 45%
```

50 time steps
```{image} ../demos/plots_rc/rc_methods_time_5_50_R_C_1.png
:width: 45%
```
```{image} ../demos/plots_rc/rc_methods_phase_5_50_R_C_1.png
:width: 45%
```

By analyzing the results, it becomes clear that the explicit Euler is underdamping, gives higher values first and then slowly decays to the stationary oscillating solution. The implicit Euler will give lower values and also decay to the stationary oscillating solution. As for the Crank Nicolson, we don't see numerical damping so we only see the stationary solution which is oscillating.
For larger time steps, for example $N=200$, we already see that the amplitudes are increasing and we don't fully recover the oscillations, which completely die out when we only use $N=50$.
So for small time steps, the Crank Nicholson recovers the behavior of the true solution the most, and our time steps should in general be small enough s.t we can recover the oscillations.

Next we look into the second case, where $RC = 10^{-4}$, therefore the homogeneus solution is now modelled by a "stiff ODE", because its transient response is vanishing really quick. This has some implications on the different methods. In general, the Explicit Euler should not be used for such models due to the nature of its stability function. For such high eigenvalues, we need really small time steps to be in a stable regime. In that case, we should expect exploding solutions for the Explicit Euler for $N \leq 10^4/2$. For the other two methods, they should remain always stable, since $\lambda \leq 0$. This time, we only simulate for $1s$. and the solution should be dominated by the cosine term this time, with the explonential decaying "quasi immediately".



5000 time steps
```{image} ../demos/plots_rc/rc_methods_time_1_10000_R_C_stiff.png
:width: 45%
```
```{image} ../demos/plots_rc/rc_methods_phase_1_10000_R_C_stiff.png
:width: 45%
```

1000 time steps with Explicit Euler
```{image} ../demos/plots_rc/rc_methods_time_1_1000_R_C_stiff.png
:width: 45%
```
```{image} ../demos/plots_rc/rc_methods_phase_1_1000_R_C_stiff.png
:width: 45%
```

1000 time steps without Explicit Euler
```{image} ../demos/plots_rc/rc_methods_time_1_1000_R_C_stiff_woEE.png
:width: 45%
```
```{image} ../demos/plots_rc/rc_methods_phase_1_1000_R_C_stiff_woEE.png
:width: 45%
```

200 time steps
```{image} ../demos/plots_rc/rc_methods_time_1_200_R_C_stiff_woEE.png
:width: 45%
```
```{image} ../demos/plots_rc/rc_methods_phase_1_200_R_C_stiff.png
:width: 45%
```

110 time steps
```{image} ../demos/plots_rc/rc_methods_time_1_110_R_C_stiff_woEE.png
:width: 45%
```
```{image} ../demos/plots_rc/rc_methods_phase_1_110_R_C_stiff.png
:width: 45%
```

So we see how for still small time steps like $\delta t = 0.001$, the Explicit Euler is already exploding and therefore should not be used for such problems.
As for the other methods they show good behavior and for larger time steps we see some unnatural oscillations in the beginning for the Crank Nicholson method, while the implicit Euler gives more accurate results since it is A and L stable. 
Again for larger time steps, the oscillations are not fully resolved and we loose accuracy.