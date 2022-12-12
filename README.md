# Longitudinal Asymmetry

In heavy-ion collisions whenever two nuclei collide, we expect an equal number of nucleons to participate from either nuclei, but it can be unequal due to nuclear density fluctuations. This gives rise to an asymmetry known as <b>longitudinal asymmetry</b>. The information of the event asymmetry allows us to isolate and study the effect of longitudinal asymmetry on various $\eta, \phi, pT$ distributions of final state particles. Here, we have looked at how fluctuations of this quantity affected various observables.

Our main aim for this project is to see how fluctuations of the parameter longitudinal asymmetry affects the distribution of various observables like pseudorapidity, rapidity and transverse momentum. In reference to https://arxiv.org/pdf/1710.07975.pdf, we have seen how they showed the dependence of asymmetry fluctuations on $\eta$ distributions for Pb-Pb collision events at $\sqrt{s_{NN}}=2.76$ TeV in the ALICE experiment conducted in 2010. For this project, we conduct simulations and generate events using the PYTHIA event generator in C++, and use the ROOT software for the analysis, plotting the relevant distributions etc; similar to what they did in https://arxiv.org/pdf/1608.01428.pdf using the HIJING MCG.

## What is Longitudinal Asymmetry?
In collisions of heavy ions, the number of nucleons participating from each of the two colliding nuclei is finite, and will fluctuate event-by-event. The geometrically overlapping region created by interacting nucleons from each nucleus is called the participant zone; and the kinematic centre of mass of the participant zone is defined as the overlap region of the colliding nuclei. In general, even at fixed
impact parameter, the number of participating nucleons
from each nucleus fluctuates around the mean due to
fluctuations in the positions of the nucleons around the
mean nuclear density profile; and event-by-event the kinematic centre of mass has a finite momentum. This finite momentum causes a <b>longitudinal asymmetry</b> in the collision; which also corresponds to a shift in rapidity, denoted as $y_0$.

When the number of nucleons participating from the two of the nuclei are $A$ and $B$, longitudinal asymmetry can be denoted as,
$
    \alpha_{part}=\frac{A-B}{A+B}
$
Here $\alpha_{part}$ denotes the asymmetry of the participant nuclei.
\subsection{Shift of rapidity}
Like we discussed before, the net finite momentum corresponds to a shift of rapidity in the participant zone, approximated as,
\begin{equation}
    y_0 \approx \frac{1}{2}ln\frac{A}{B}
\end{equation}

\pagebreak
w.r.t how we defined $\alpha_{part}$ in \ref{longasymm}, we can rewrite $y_0$ in terms of $\alpha_{part}$ as,
\begin{equation}
    y_0=\frac{1}{2}ln\frac{1+\alpha_{part}}{1-\alpha_{part}}
\end{equation}
The unequal number of nucleons in the participant zone implies an unequal number of spectators of the two colliding nuclei, $N-A$ and $N-B$ , respectively, where N is the total number of nucleons in each nucleus. The spectator asymmetry,
\begin{equation}
    \alpha_{spec}=\frac{(N-A)-(N-B)}{(N-A)+(N-B)}
\end{equation}
is related to the participant asymmetry by
\begin{equation}
    \alpha_{spec}=-\alpha_{part}\frac{A+B}{2N-(A+B)}
\end{equation}
and $y_0$ is related to the spectator asymmetry by,
\begin{equation}
    y_0=\frac{1}{2}ln\frac{(A+B)(1+\alpha_{spec})-2N\alpha_{spec}}{(A+B)(1+\alpha_{spec})+2N\alpha_{spec}}
\end{equation}
In our present work, we have generated 10,000 events of $0-20\%$ and $80-100\%$ centralities of Pb-Pb collisions at $\sqrt{s_{NN}}=2.76$ TeV using PYTHIA.

