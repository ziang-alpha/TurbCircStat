rsConsider the integration:
$$
\langle\varGamma_C^2\rangle_\text{eq} = \iint \frac{1}{\alpha + \beta(k_1^2+k_2^2)}\frac{\sin^2(k_1l_1)\sin^2(k_2l_2)}{k_1^2k_2^2}\mathrm{d}k_1\mathrm{d}k_2
$$
Apply the Laplacian transformation:
$$
\begin{aligned}
\langle\varGamma_C^2\rangle_\text{eq} =& \iint \int_0^\infty e^{-t(\alpha + \beta(k_1^2+k_2^2))}\mathrm{d}t\frac{\sin^2(k_1l_1)\sin^2(k_2l_2)}{k_1^2k_2^2}\mathrm{d}k_1\mathrm{d}k_2\\
=& \int_0^\infty e^{-t\alpha} \left(\int \frac{\sin^2(k_1l_1)}{k_1^2 e^{t\beta k_1^2}}\mathrm{d}k_1\right)\left(\int \frac{\sin^2(k_2l_2)}{k_2^2 e^{t\beta k_2^2}}\mathrm{d}k_2\right)\mathrm{d}t
\end{aligned}
$$

Here
$$
\begin{aligned}
\int \frac{\sin^2(k_1l_1)}{k_1^2 e^{t\beta k_1^2}}\mathrm{d}k_1 = & l_1\int \frac{\sin^2(k_1l_1)}{(k_1l_1)^2 }e^{-t\beta l_1^{-2} (k_1l_1)^2}\mathrm{d}k_1l_1\\ 
    =& l_1 \int \frac{\sin^2x}{x^2}e^{-t\beta l_1^{-2}x^2}\mathrm{d}x
\end{aligned}
$$

The integration over the whole $\mathbb{R}$ reads:
$$
l_1 \int_{-\infty}^\infty \frac{\sin^2x}{x^2}e^{-t\beta l_1^{-2}x^2}\mathrm{d}x = l_1\left(\pi\operatorname{erf}\left(\frac{l_1}{\sqrt{\beta t}}\right)-\frac{\sqrt{\pi\beta t}}{l_1}(1-\exp(-l_1^2/\beta t))\right):=l_1F(\sqrt{\beta t}/l_1)
$$

So
$$
\begin{aligned}
\langle\varGamma_C^2\rangle_\text{eq} =& l_1l_2\int_{0}^\infty e^{-\alpha t}F(\sqrt{\beta t}/l_1)F(\sqrt{\beta t}/l_2)\mathrm{d}t\\
    =&A\alpha^{-1}\int_{0}^\infty e^{-\tau}F(l_\text{eq}\sqrt{\tau}/l_1)F(l_\text{eq}\sqrt{\tau}/l_2)\mathrm{d}\tau
\end{aligned}
$$