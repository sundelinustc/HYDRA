# HYDRA
We applied a novel non-linear learning algorithm termed HYDRA (Heterogeneity through Discriminative Analysis) developed by Dr. Sotiras and colleagues to do binary classification and subtype identification simultaneously. HYDRA is different from the widely used supervised learning techniques including support vector machines and random forests that are unable to distinguish between subtypes of cases without training, and is also different from the popular unsupervised clustering techniques including k-means and community detection that are vulnerable to confounding inter-individual variations that are irrelevant to disease (e.g. age, sex). HYDRA clusters patients based on their differences from controls, and is thus more appropriate to discover phenotypic heterogeneity of underlying neurobiological processes. 
![alt text](https://github.com/sundelinustc/HYDRA/blob/Figure1.png?raw=true)

To accomplish this, HYDRA finds multiple hyperplanes, which together form a convex polytope (polygon in n-dimensional space) that separates controls from subtypes. Rather than coercing participant data points into a single common discriminative pattern, HYDRA allows for the separation of groups distinguished by multiple decision boundaries. The result is a data-driven approach to identifying disease subtypes that can be evaluated further with independently-measured clinical and imaging characteristics. HYDRA allows covariates to be introduced during the clustering procedure. This approach is potentially useful because it seeks to identify subtypes within cases, which are independent of symptoms, but more critically the subtypes are neurobiologically distinct. 
![alt text](https://github.com/sundelinustc/HYDRA/blob/Figure3.png?raw=true)

In a summary, HYDRA is a hybrid between unsupervised clustering and supervised classification methods by simultaneously fitting maximum margin classification boundaries and elucidating disease subtypes, which has been validated on simulated and clinical data.
