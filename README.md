# catalytic_model_US


This code reproduces the analysis in the paper ``Age-specific case data reveal varying dengue transmission intensity in US states and territories''. The code runs with simulated population census and case data to run the catalytic model.

Three models exist:
- Model S: 

$$I_2(a, t) =  3  \lambda (t) \text{ . } M(a, t), $$

For the **S** model, the age-specific incidence of secondary infections, $I_2$, for an individual who already experienced a monotypic infection is then the force of infection of the 3 other serotypes that individual can be infected with:


- Model PS {-}
The model described above was based on cases reported being exclusively secondary infections (model **S**). For the second model, **PS**, we assumed that reported cases could be both primary and secondary infections, the equation \@ref(eq:incidence) for the age-specific incidence is now:


$$ I_{12}(a, t) = 4  \lambda(t) \text{ . } S(a, t) + 3 \lambda(t) \text{ . } M(a, t)$$

- Model P {-}
The last model assumes that reported cases are exclusively primary infections. The age-specific incidence at time $t$ is then:

$$ I_1(a, t) =  4  \lambda(t) \text{ . } S(a, t). $$

