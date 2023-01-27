# catalytic_model_US


This code reproduces the analysis in the paper "Age-specific case data reveal varying dengue transmission intensity in US states and territories''. The code runs a catalytic model with simulated population census and case data.

Use ``code/1-main.R`` to run the catalytic models. Priors and saved outputs can be modified in the ``code/0-export_code.R``.

We explored three different models, with corresponding stan code called. Models are named "model S", "model P" or "model PS" following:

- **Model S:** 
The age-specific incidence of secondary infections, $I_2$, for an individual who already experienced a monotypic infection is the force of infection of the 3 other serotypes that individual can be infected with:

$$I_2(a, t) =  3  \lambda (t) \text{ . } M(a, t), $$

- **Model PS:**
For model **PS**, it is assumed that reported cases could be both primary and secondary infections, age-specific incidence at time $t$ is now:


$$ I_{12}(a, t) = 4  \lambda(t) \text{ . } S(a, t) + 3 \lambda(t) \text{ . } M(a, t)$$

- **Model P**
This last model assumes that reported cases are exclusively primary infections. The age-specific incidence at time $t$ is then:

$$ I_1(a, t) =  4  \lambda(t) \text{ . } S(a, t). $$
