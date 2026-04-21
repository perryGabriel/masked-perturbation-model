# Future-work roadmap for the masked-perturbation program

This document consolidates the open questions from:

- **Chapter 3, Section V** of the main masked-perturbation manuscript
- **Chapter 4, Future Work and Conclusion** of the future-work PDF
- the practical limitations called out at the end of Chapter 4
- one additional **nonlinear ISS extension** that we discussed separately because it fits the same research line

The goal is not just to list questions, but to rank them, explain which are promising or weak, sketch likely solution paths, and provide **Codex-ready implementation prompts** for the parts that are computational.

---

## 0. Executive prioritization

### Tier 1: most promising, high-leverage, tractable enough for near-term progress
1. **Simultaneous destabilization / dominant attacks across multiple defense masks**
2. **Mixed-defense action-set design via success-set geometry**
3. **Stackelberg / leader-follower version of the masked game**
4. **Signal-importance ranking and approximation algorithms for defense construction**
5. **Expected-vulnerability reformulation of the full FMP under random exposure**

### Tier 2: promising, but probably second-wave work
6. **Alternative defense objectives and trade-off frontiers**
7. **When it is not optimal to spend the whole budget**
8. **Defense design that intentionally forces attacker mixing**
9. **Sensitivity of the optimal defense to the attack action**
10. **Single-defense coverage / set-cover style approximation guarantees**

### Tier 3: important, but likely theory-heavy or too easy to become ill-posed
11. **Existence/minimality certificates for simultaneous attacks using interpolation**
12. **Partial/probabilistic/delayed masking as uncertainty sets**
13. **Time-varying schedules, saturations, packet drops, and non-LTI channels**
14. **The original “hide or rotate the dynamics away from exposed signals” full-FMP formulation**

### Tier 4: adjacent but strategically valuable
15. **Nonlinear ISS / small-gain extension of masked perturbations**

---

## 1. Source map: where each question came from

### From Chapter 3, Section V of the main text
- **V-A. Simultaneous Destabilization (Attacker’s Perspective)**:
  - minimal one-attack destabilization of multiple masked plants
  - simultaneous-destabilization conditions
  - remarks about subset/superset relationships and shared footprints
- **V-B. Stackelberg Security Game (Defender’s Perspective)**:
  - defender commits first
  - attacker observes and responds
  - resource-allocation interpretation

### From Chapter 4 of the future-work PDF
#### Attacker-side questions
- whether minimal destabilizing attacks for different masks simultaneously destabilize multiple masks
- necessary/sufficient dominance conditions from subset/superset structure or common vulnerable links
- existence of dominant admissible attacks for any set of `k` defenses
- whether dominant attacks are minimal
- Nevanlinna–Pick-style interpolation constructions
- impossibility certificates
- mixed-attack action-set construction via success-set overlap / set cover / intersection-over-union
- approximation-ratio guarantees

#### Defender-side questions
(i) better optimization formulations for defense policy  
(ii) sensitivity of optimal defense to the attack  
(iii) sub-exponential algorithm for constructing one defense that covers an attack  
(iv) approximation guarantees if exact construction is intractable  
(v) whether it is always optimal to spend the full budget  
(vi) whether the defender can force attacker mixing  
(vii) mixed defense construction when no single defense covers all attacks  
(viii) complexity / approximation guarantees for that mixed-defense construction  
(ix) sub-exponential methods for determining which signals are most important to secure

#### Reconsidering the full FMP
- minimizing `Grw` by realization design
- why “hiding/rotating” dynamics away from known exposed channels is probably not a satisfying formulation
- replacing that with **expected vulnerability over random exposed subsets**

#### Limitations explicitly named in Chapter 4
- binary masks idealize coverage; extension to partial/probabilistic/delayed enforcement
- simultaneous stabilization/destabilization remains hard
- LTI/Laplace setting excludes time variation, saturation, packet drops, and queuing

---

# 2. Detailed ranked plan

---

## 2.1 Simultaneous destabilization / dominant attacks across multiple masks
**Priority:** 1  
**Assessment:** very promising  
**Why it is good:** This is the most natural next attacker-side problem and it is already partly scaffolded by the remarks in Chapter 3. The masked plants differ only by binary row/column suppression, so this is more structured than generic simultaneous stabilization. That gives you leverage.  
**Main source questions covered:** Ch. 3 Sec. V-A; Ch. 4 attacker paragraph 1.

### Core question
Given defense masks `∇_1, ..., ∇_k`, can one admissible attack `Δ` destabilize all corresponding masked maps `M_{∇_j}` simultaneously? If yes:
- when does such a `Δ` exist,
- when is it minimal,
- how large must it be,
- when can we certify nonexistence?

### Why it is likely publishable
It sits directly between:
- robust-control simultaneous stabilization,
- structured attack-surface masking,
- game-theoretic dominance of attack actions.

The Chapter 3 remarks already suggest exploitable structure:
- if `∇_i ⊆ ∇_j`, then an attack successful against `∇_j` may remain successful against `∇_i`
- shared/unshared accessible footprints suggest necessary combinatorial conditions
- disjoint footprints admit direct-sum constructions

### Likely theorem targets
1. **Necessary support-overlap condition**  
   If two defense profiles expose no common read/write channels, then any simultaneous attack must decompose additively across the two footprints.

2. **Sufficient condition from subset ordering**  
   If the defended sets are nested, then a minimal attack for the most restrictive mask may dominate weaker masks under admissibility assumptions.

3. **Existence criterion under unique vulnerable frequencies**  
   If each masked plant has a destabilizing pressure point at distinct frequencies and the admissible attack class contains a stable interpolant meeting those constraints, then a simultaneous destabilizer exists.

4. **Impossibility certificate**  
   No admissible `Δ` exists if the interpolation constraints conflict with admissibility (stable/proper/order/budget) or if footprint support constraints are incompatible.

### Attempted solution path
#### Theory
- Start with `k=2`, scalar or diagonal-structured attacks.
- Use the Chapter 3 remarks as lemmas.
- Re-express simultaneous destabilization as either:
  - a structured interpolation problem, or
  - a feasibility problem over frequency-sampled pressure constraints.
- Prove simple support-based necessary conditions before trying minimality.

#### Computation
- For each mask:
  1. compute a candidate minimal destabilizing `Δ_j`
  2. test pairwise cross-success on other masks
  3. search for additive or interpolated combinations
- Build a graph whose nodes are masks and whose edges indicate shared vulnerable footprint or shared vulnerable frequency.
- Use this graph to identify candidate groups that may admit simultaneous attacks.

### Risks
- Exact minimality may be very hard beyond `k=2`.
- Simultaneous destabilization may become a hard interpolation problem very quickly.
- “Dominant” may depend strongly on the admissible attack class.

### Recommended background to look into
- Zhou, Doyle, Glover — *Robust and Optimal Control* (small gain, LFT)
- Blondel — *Simultaneous stabilization of linear systems*
- Cui & Lindquist (2023, 2024) on simultaneous stabilization via analytic interpolation
- Doyle, Francis, Tannenbaum — *Feedback Control Theory*
- Any literature on structured singular value with multiple plants

### Concrete outputs
- a theorem/proposition section for `k=2`
- a numerical catalog of pairwise simultaneous destabilization outcomes
- a small benchmark dataset of masked plants and candidate simultaneous attacks
- a reduced-game attack-action generator based on simultaneous attacks

### Codex prompt for code work
```text
Implement a research prototype for simultaneous destabilization in the masked-perturbation package.

Scope:
1. Accept a base model map M and a finite set of defense masks ∇_1, ..., ∇_k.
2. For each mask, compute or accept a candidate minimal destabilizing attack Δ_j.
3. Build utilities to test whether a given Δ destabilizes each masked plant M_{∇_j}.
4. Compute:
   - pairwise cross-success matrix
   - support-overlap matrix
   - candidate dominance relations among attacks
   - additive combined attacks Δ_i + Δ_j
5. Add a frequency-sampled interpolation-based search:
   - identify “pressure frequencies” for each mask
   - fit a stable proper scalar or diagonal transfer function Δ that approximately matches prescribed gains/phases at those frequencies
   - test simultaneous destabilization success afterward
6. Produce reports:
   - masks covered by each attack
   - attack norms
   - simultaneous-attack success tables
   - figures showing vulnerable frequencies and interpolation fits

Important:
- be honest that this is an exploratory search, not a complete simultaneous-destabilization solver
- support at least scalar, diagonal, and small 2x2 examples
- keep the public API small and testable
```

### Completion update (implemented)
- Added module: `src/mpmgame/simultaneous_attacks.py` with data models (`MaskDefense`, `AttackCandidate`, `SimultaneousProblem`), per-mask success evaluation, cross-success and support-overlap matrices, dominance detection, additive-combination search, and output writers.
- Exposed the simultaneous-destabilization API from `src/mpmgame/__init__.py`.
- Added runnable demo: `scripts/simultaneous_destabilization_demo.py`.
- Added notebook: `notebooks/simultaneous_destabilization_demo.ipynb`.
- Added report: `reports/simultaneous_destabilization/report.md`.
- Added tests: `tests/test_simultaneous_destabilization.py`.

Implemented command examples:
```bash
python scripts/simultaneous_destabilization_demo.py
pytest tests/test_simultaneous_destabilization.py
```

Produced artifacts under `results/simultaneous_destabilization/`:
- `cross_success_matrix.csv`
- `support_overlap_matrix.csv`
- `dominance_relations.json`
- `coverage_by_attack.json`
- `additive_attacks.json`
- `run_summary.json`
- `cross_success_heatmap.png`
- `support_overlap_heatmap.png`
- `coverage_frequency.png`

---

## 2.2 Mixed-defense action-set design via success-set geometry
**Priority:** 2  
**Assessment:** very promising  
**Why it is good:** This is one of the cleanest algorithmic questions in the document because the success-set formulation is already explicit. It naturally connects to maximum coverage, set cover, submodularity, and approximation algorithms.  
**Main source questions covered:** defender (vii), (viii), (ix); also practical-takeaway “use success sets to prune”.

### Core question
When no single defense `∇` can stabilize against all admissible attacks, can we choose a **small bounded subset** of defenses whose mixed strategy performs nearly as well as the full exponential defense space?

The future-work document even proposes a surrogate objective:

$$
\zeta^* = \arg\max_{\zeta \subseteq \nabla}
-\alpha_z |\zeta|
+\sum_{i,j}\alpha_0 |S(\nabla_i)\cup S(\nabla_j)|
+\alpha_1 |S(\nabla_i)\cap S(\nabla_j)|
$$

subject to `|ζ| ≤ z` and budget feasibility.

### Why it is strong
- computationally concrete
- immediate path to experiments
- likely admits approximation guarantees if the objective is shown submodular or approximately submodular
- useful even if the deeper simultaneous-destabilization theory is incomplete

### Attempted solution path
#### Theory
1. Formalize each defense as a set `S(∇) ⊆ Δ` of attacks it defeats.
2. Define a coverage objective:
   - expected success under a prior on attacks, or
   - worst-case support coverage, or
   - the proposed union/intersection-weighted objective.
3. Check whether the objective is monotone submodular or reducible to known submodular families.
4. If yes, greedy gives a `(1 - 1/e)`-style guarantee for the set-selection stage.
5. Then solve the reduced zero-sum game over the selected defense set.

#### Computation
- Enumerate or sample admissible defenses up to a budget.
- Compute success sets.
- Run:
  - greedy maximum coverage
  - greedy submodular selection
  - beam search
  - random-restart local search
- Compare reduced-set mixed strategies with the full game when tractable.

### Risks
- the exact union/intersection objective may not be submodular
- good reduced sets for `u_a = 1, u_d = 0` may differ from good reduced sets under richer utilities
- pairwise overlap metrics can overvalue redundant defenses

### Recommended background
- Nemhauser-Wolsey style submodular maximization
- maximum coverage / budgeted maximum coverage
- Stackelberg security games with combinatorial action spaces
- set cover and hitting set approximations
- the security-game references already listed in Chapter 4

### Concrete outputs
- approximation algorithm paper/section
- benchmark comparing full-game vs reduced-action-set performance
- empirical curves: value vs action-set size `z`
- ablation on the `α_z, α_0, α_1` weights

### Codex prompt
```text
Extend the masked-perturbation package with reduced mixed-defense action-set design.

Implement:
1. Success-set computation S(∇) over a finite attack set.
2. Defense-set selection algorithms under a cardinality cap z:
   - greedy maximum coverage
   - greedy union/intersection objective
   - beam search
   - local search
3. Optional attack prior p(Δ), and support expected-coverage objectives.
4. After choosing a defense subset ζ, solve the reduced zero-sum game on ζ.
5. Compare reduced-game value against the full game whenever the full game is small enough.

Deliver:
- plots of value vs z
- coverage-vs-z plots
- overlap heatmaps of success sets
- markdown report summarizing how close the reduced defense subset comes to full-game performance

Important:
- label any approximation guarantee only when it is actually justified by the objective used
- keep code modular so alternative objectives can be swapped in
```

### Completed implementation note
- Added mixed-defense module: `src/mpmgame/mixed_defense.py` (success-set extraction, capped subset selection with greedy/random baselines, and union/intersection/cardinality objective components).
- Added reduced-game subset evaluation helpers in `src/mpmgame/core.py` and exported the APIs in `src/mpmgame/__init__.py`.
- Added benchmark sweep runner: `scripts/mixed_defense_benchmark.py` with artifacts saved to `results/mixed_defense/` and `reports/mixed_defense/`.
- Added notebook: `notebooks/mixed_defense_action_set_demo.ipynb`.
- Added report: `reports/mixed_defense/report.md`.
- Added tests: `tests/test_mixed_defense.py`.

---

## 2.3 Stackelberg / leader-follower masked security game
**Priority:** 3  
**Assessment:** very promising  
**Why it is good:** The model already looks like a security game, and the Stackelberg version is operationally closer to actual deployments where the defender commits to a patrol/monitoring policy and the attacker observes it. This is a natural next paper once the normal-form mixed game is implemented.  
**Main source questions covered:** Ch. 3 Sec. V-B; Chapter 4 references on security games.

### Core question
What changes when the defender commits first to a mixed strategy over masks, and the attacker best-responds after observing that policy?

### Why it is good
- there is mature algorithmic literature to borrow from
- the attack/defense semantics are interpretable
- it may generate cleaner prescriptions than the simultaneous-move zero-sum model
- it fits the application language of guards, scans, and random monitoring

### Attempted solution path
#### Theory
- Start with exact commitment against a finite attack set.
- Compare:
  - simultaneous-move minimax value
  - defender Stackelberg value
- Analyze whether dominated-action pruning still preserves optimal commitment solutions.

#### Computation
- For finite attack/defense sets, solve the defender commitment LP / MILP.
- Support attacker tie-breaking variants:
  - strong Stackelberg equilibrium (defender-favorable tie-break)
  - weak / adversarial tie-break
- Compare committed distributions to simultaneous minimax distributions.

### Risks
- if utilities remain strictly zero-sum and finite, Stackelberg may collapse to a similar value in some formulations
- the interesting part may require observational uncertainty or bounded rationality

### Recommended background
- Tambe / Kiekintveld line of Stackelberg security game work
- strong Stackelberg equilibrium formulations
- online learning / commitment without regret in security games
- randomized resource allocation literature

### Concrete outputs
- exact Stackelberg solver for reduced games
- comparison notebook: minimax vs Stackelberg
- study of attack observability / tie-breaking assumptions

### Codex prompt
```text
Add a Stackelberg-security-game module to the masked-perturbation package.

Requirements:
1. Input: finite attack set Δ, defense set ∇, payoff matrix U.
2. Implement defender-commitment optimization for:
   - strong Stackelberg equilibrium
   - adversarial tie-breaking
3. Output:
   - defender committed mixed strategy
   - attacker best response(s)
   - defender value
4. Compare against the simultaneous-move minimax solution on the same game.
5. Add figures:
   - support of the optimal defense strategy
   - comparison of values across several example games
6. Integrate with dominated-strategy elimination and reduced games.

Be explicit in docs about tie-breaking assumptions and when the formulation is zero-sum versus general-sum.
```

### Completion update (implemented)
- Added Stackelberg module: `src/mpmgame/stackelberg.py` with defender commitment over budget-feasible masks, attacker best-response utilities for pure/mixed commitments, mixed/pure commitment solvers, and JSON/CSV export helpers.
- Exported Stackelberg APIs from `src/mpmgame/__init__.py`.
- Added runnable demo script: `scripts/stackelberg_demo.py`.
- Added notebook reproduction: `notebooks/stackelberg_masked_game_demo.ipynb`.
- Added usage guide: `docs/stackelberg.md`.
- Added consistency tests: `tests/test_stackelberg.py`.

Produced artifacts:
- `results/stackelberg/mixed_summary.json`
- `results/stackelberg/mixed_policy_table.csv`
- `results/stackelberg/pure_summary.json`
- `results/stackelberg/pure_policy_table.csv`
- `results/stackelberg/mixed_expected_attack_utilities.csv`
- `results/stackelberg/toy_payoff_matrix.csv`
- `reports/stackelberg/response_map.png`
- `reports/stackelberg/value_comparison.png`

Reproducibility commands:
```bash
python scripts/stackelberg_demo.py
pytest tests/test_stackelberg.py
```

---

## 2.4 Signal-importance ranking and approximation algorithms for defense construction
**Priority:** 4  
**Assessment:** very promising  
**Why it is good:** This is practically important and can likely be turned into strong heuristics even if exact sub-exponential algorithms are impossible.  
**Main source questions covered:** defender (iii), (iv), (ix).

### Core question
Can we rank exposed signals or small signal groups by their marginal impact on vulnerability or success sets, and use that ranking to build near-optimal defenses?

### Attempted solution path
- Define importance scores such as:
  - decrease in worst-case `μ(M_∇)` when a signal is secured
  - decrease in attack-success coverage
  - Shapley-style marginal importance over signal subsets
  - frequency-weighted score using row/column importance in `M(jω)`
- Compare:
  - greedy single-signal selection
  - pairwise lookahead
  - Shapley approximation by Monte Carlo
- Test whether these heuristics produce strong defense subsets quickly.

### Risks
- importance may be highly non-additive
- pairwise synergies can dominate individual rankings
- ranking could depend heavily on the chosen threat model

### Recommended background
- submodular feature selection
- Shapley-value attribution for combinatorial design
- robust sensor/actuator placement
- network interdiction / defender resource allocation

### Concrete outputs
- signal-importance ranking metric catalog
- benchmark of greedy rankings vs exact optimum on small problems
- practical heuristic for large spaces

### Codex prompt
```text
Implement signal-importance ranking heuristics for masked defenses.

Add:
1. Marginal vulnerability-drop score for securing each read/write signal.
2. Success-set marginal gain score under a finite attack set.
3. Pairwise and optional Monte Carlo Shapley-style scores.
4. Greedy defense construction using these scores under a cost budget.
5. Benchmark exact-vs-greedy performance on small games where exhaustive search is possible.

Deliver:
- ranked signal tables
- plots of cumulative value vs budget
- exact-vs-greedy gap plots
- support for threat-model-specific scoring
```

---

## 2.5 Expected-vulnerability reformulation of the full FMP under random exposure
**Priority:** 5  
**Assessment:** promising, but only under the **reformulated** version  
**Why it is good:** The document itself says the naive formulation “hide or rotate the dynamics away from known exposed signals” is not very promising. The random-exposure expected-vulnerability reformulation is the right salvage.  
**Main source questions covered:** 4.1.3 Reconsidering the Full FMP.

### Split assessment
#### Bad version
**Question:** choose a realization so that `M = G_rw` is tiny at known exposed channels.  
**Why it is weak:** if the access points are known and independent of the contract, the defender can artificially hide vulnerable dynamics in unexposed coordinates. That makes the problem brittle, access-model-dependent, and somewhat “honeypot-like.”

#### Good version
**Question:** choose a realization to minimize **expected vulnerability over a distribution of possible exposed subsets**, with every state/signal having nonzero exposure probability.  
**Why it is better:** it blocks the trivial hiding solution and turns the problem into a meaningful stochastic design problem.

### Attempted solution path
#### Theory
- Parameterize admissible realizations of a fixed behavioral contract.
- Define a distribution `π(E)` over exposed read/write sets.
- Objective:

$$
\min_{\text{realization}} \mathbb{E}_{E \sim \pi}\left[V(M_E)\right]
$$

  or a risk-sensitive version like CVaR.
- Study whether separability survives in expectation.

#### Computation
- Restrict to small DSF/state-space parameterizations.
- Sample exposure sets from `π`.
- Use stochastic search / projected gradient / surrogate optimization over realization parameters.
- Compare:
  - expected vulnerability
  - worst-case vulnerability
  - variance / CVaR of vulnerability

### Risks
- realization parameterization may be awkward
- structural identifiability constraints may dominate the optimization
- exact constraints for preserving `G_yu` may be difficult

### Recommended background
- DSF / structured realization literature
- stochastic robust design
- risk-sensitive design (CVaR, distributionally robust optimization)
- sensor/actuator exposure models
- the FMP appendix and DSF references already in your thesis sources

### Concrete outputs
- a precise problem statement
- a toy stochastic-realization experiment
- comparison between “known fixed exposure” and “random exposure” objectives

### Codex prompt
```text
Implement a toy research module for the full-FMP random-exposure reformulation.

Scope:
1. Work only on small systems with a fixed nominal behavioral map G_yu.
2. Parameterize a family of realizations (for example, DSF-style or small state-space realizations that preserve the same external contract).
3. Define a distribution over exposed read/write subsets, with every signal having nonzero exposure probability.
4. For each realization:
   - sample many exposure subsets
   - construct the corresponding exposed map M_E
   - compute vulnerability statistics: mean, worst-case, variance, optional CVaR
5. Implement a simple search over realization parameters:
   - random search
   - coordinate search
   - optional evolutionary or Bayesian optimization
6. Produce:
   - distributions of vulnerability under exposure uncertainty
   - Pareto plots of expected vs worst-case vulnerability
   - examples where fixed known exposure gives misleadingly optimistic designs

Important:
- be explicit that this is a toy expected-vulnerability study, not a full solver for the original FMP
```

---

## 2.6 Alternative defense objectives and trade-off frontiers
**Priority:** 6  
**Assessment:** good  
**Main source questions covered:** defender (i)

### Core question
Is Eq. (4) the right optimization problem, or should the defender instead:
- minimize cost subject to risk,
- minimize CVaR of attack success,
- minimize worst-case `μ`,
- maximize expected stabilization probability under uncertainty?

### Attempted solution path
- Define a family of defender problems:
  1. fixed-budget risk minimization
  2. minimum-cost for acceptable risk
  3. CVaR / entropic risk
  4. robust-to-uncertainty-set variants
- Compare resulting strategies on the same reduced games.

### Why it is good
This is easy to turn into a clean computational section and may clarify which objective is most operational.

### Risks
Mostly empirical unless tied to a real decision-theoretic argument.

### Recommended background
- risk-sensitive optimization
- CVaR in security games
- robust vs stochastic optimization trade-offs

### Codex prompt
```text
Add alternative defender objectives to the masked game package.

Implement and compare:
1. minimize expected attacker utility under fixed budget
2. minimize defense cost subject to expected attacker utility <= epsilon
3. minimize CVaR of attacker success
4. optional max-min robust version over uncertainty in attack probabilities

Generate Pareto frontiers:
- cost vs expected destabilization probability
- cost vs CVaR
- budget vs game value

Return defense supports and value tables for each objective.
```

---

## 2.7 When is it not optimal to spend the whole budget?
**Priority:** 7  
**Assessment:** good and surprisingly important  
**Main source questions covered:** defender (v)

### Core question
Must an optimal defender always exhaust its budget?

### Why it matters
If the answer is “yes under broad conditions,” then you can prune the action space to Pareto-front profiles. If “no,” many current simplifications become invalid.

### Attempted solution path
- Construct small counterexamples first.
- Look for sufficient conditions under which spending more never hurts:
  - monotone attack-success structure
  - monotone threat models
  - utilities with no footprint-narrowing paradox
- Then characterize failure modes:
  - narrowing footprint makes the attacker’s best response more concentrated
  - mixed-strategy geometry can make “over-defending” bad

### Risks
This may produce many edge-case counterexamples and few elegant theorems.

### Recommended background
- monotone comparative statics
- resource-allocation paradoxes in security games
- congestion/coverage paradoxes

### Concrete outputs
- theorem under monotonicity assumptions
- counterexample notebook outside those assumptions

### Codex prompt
```text
Build a search tool for budget-monotonicity counterexamples in the masked-defense game.

Requirements:
1. Enumerate finite attack/defense games generated from small masked systems.
2. For each budget level C, compute the defender optimum.
3. Detect whether increasing budget strictly worsens defender value or changes support in non-monotone ways.
4. Save any counterexample instances with:
   - payoff matrix
   - action labels
   - optimal strategies by budget
   - explanation plots

Also test sufficient monotonicity conditions on restricted game families.
```

---

## 2.8 Can the defender force the attacker to mix?
**Priority:** 8  
**Assessment:** good  
**Main source questions covered:** defender (vi)

### Core question
Can the defender deliberately choose a mask or mask distribution that destroys dominant attacks and forces the attacker to spread probability mass, thereby increasing the chance of stable pairings?

### Attempted solution path
- Formalize “force mixing” as:
  - increase in attacker support size
  - decrease in maximal pure-action success probability
  - creation of multiple equal-value attacker best responses
- Search for design principles:
  - symmetry
  - balancing vulnerable frequencies
  - eliminating dominant success-set inclusions

### Risks
Depends strongly on attacker observability and tie-breaking.

### Recommended background
- entropy-regularized games
- security games with deceptive signaling / commitment
- support-size and sparsity in equilibrium

### Codex prompt
```text
Add a “force attacker mixing” analysis module.

Implement:
1. metrics of attacker support size and concentration
2. search over defense actions or defense subsets that maximize these metrics subject to budget
3. comparison between:
   - defender value
   - attacker support entropy
   - stable-pairing probability
4. visualize how defense design changes the attacker best-response correspondence
```

---

## 2.9 Sensitivity of the optimal defense to the attack
**Priority:** 9  
**Assessment:** moderate but useful  
**Main source questions covered:** defender (ii)

### Core question
For a single attack `Δ`, how does the optimal defense `∇*` change with `Δ`? Is there a meaningful analogue of `∂∇*/∂Δ`?

### Assessment
Good diagnostic question, but probably not the highest-priority paper because the defense is discrete. The right object is likely **piecewise-constant policy regions**, not a derivative in the classical sense.

### Attempted solution path
- Treat `Δ` as parameterized by a low-dimensional vector (gain, pole, support).
- Partition parameter space into regions where the same optimal defense remains optimal.
- Study boundary crossings and tie regions.

### Risks
Derivative notation is misleading because `∇*` is discrete.

### Recommended background
- parametric combinatorial optimization
- piecewise-constant policy maps
- best-response region analysis

### Codex prompt
```text
Implement a parameter-sweep study of defense sensitivity to a parameterized attack Δ(θ).

Requirements:
1. define low-dimensional attack families Δ(θ)
2. for each θ on a grid, solve for optimal defense
3. identify regions of constant optimal defense
4. plot:
   - θ vs optimal defense label
   - θ vs defender value
   - transitions / tie points
5. support scalar and small MIMO examples
```

---

## 2.10 Single-defense coverage and approximation guarantees
**Priority:** 10  
**Assessment:** moderate-good  
**Main source questions covered:** defender (iii), (iv)

### Core question
Can we construct a single defense `∇` that covers a target attack or a large expected-success set efficiently?

### Assessment
Good as a foundational subproblem, but it may collapse into a version of weighted set cover / hitting set. That is still useful, but perhaps not by itself a full paper unless the structure of masks yields something special.

### Attempted solution path
- Formalize the exact decision problem.
- Determine whether it reduces to:
  - set cover,
  - hitting set,
  - maximum coverage,
  - knapsack with set-valued rewards.
- Then derive:
  - NP-hardness by reduction if appropriate
  - greedy approximation where possible
  - exact branch-and-bound for small cases

### Recommended background
- set cover / hitting set
- budgeted maximum coverage
- combinatorial optimization on bipartite incidence graphs

### Codex prompt
```text
Implement exact and approximate solvers for single-defense coverage.

Input:
- finite attack set
- admissible defense set with costs
- success sets S(∇)

Tasks:
1. exact search for the cheapest defense covering a target subset of attacks
2. greedy approximation for maximum covered attack mass under budget
3. branch-and-bound for small instances
4. output approximation gaps whenever exact optimum is known
```

---

## 2.11 Interpolation conditions, minimality, and impossibility certificates
**Priority:** 11  
**Assessment:** promising but theory-heavy and risky  
**Main source questions covered:** attacker interpolation paragraph

### Core question
Can Nevanlinna–Pick-style interpolation characterize existence of a simultaneous destabilizing `Δ`, and can it certify impossibility?

### Assessment
This is mathematically deep and potentially elegant, but it is probably **not** the best next project unless you already have the right interpolation machinery lined up. It is a good second-paper theory direction once the simultaneous-destabilization numerics are mature.

### Attempted solution path
- Reduce simultaneous destabilization to interpolation constraints at selected frequencies.
- Start with scalar stable proper `Δ`.
- Determine whether the needed interpolation conditions fit classical Schur / Nevanlinna–Pick formulations or require derivative / weighted constraints.
- Use failure of the Pick matrix or analogous feasibility object as an impossibility certificate.

### Risks
- easy to sink a lot of time into function-theoretic details
- may require stronger regularity assumptions than the applied problem naturally provides

### Recommended background
- classical Nevanlinna–Pick interpolation
- simultaneous stabilization via analytic interpolation (Cui & Lindquist)
- strong stabilization and rational interpolation

### Code note
This is not a primary Codex target. Code should support experimentation only, not theorem discovery.

---

## 2.12 Partial/probabilistic/delayed masking as uncertainty sets
**Priority:** 12  
**Assessment:** good long-term extension  
**Main source questions covered:** limitation 1 in Chapter 4

### Core question
How should binary masks be generalized when enforcement is partial, probabilistic, delayed, or uncertain?

### Assessment
Good and practically important, but it opens a much broader robust-control problem. I would postpone it until the binary-mask program is mature.

### Attempted solution path
- Replace binary mask entries by:
  - probabilities,
  - attenuation factors,
  - uncertain bounded parameters,
  - switching processes.
- Define corresponding uncertainty set for the masked map or for the implemented coverage policy.
- Compare:
  - mixed strategies over binary masks
  - true uncertainty-set robust design

### Risks
- can easily blur the clean difference between randomized policies and parametric uncertainty
- may require stochastic or hybrid systems tools

### Recommended background
- uncertain systems / μ-analysis
- Markov jump linear systems
- robust switching / packet-drop models
- security games under execution uncertainty

### Code note
A small Monte Carlo module would be fine, but avoid overcommitting to this before the core binary-mask theory is settled.

---

## 2.13 Time-varying schedules, saturations, packet drops, queues, and non-LTI channels
**Priority:** 13  
**Assessment:** important but broad  
**Main source questions covered:** limitation 3 in Chapter 4

### Assessment
This is a broad family of extensions rather than a single crisp question. Good eventual direction; weak immediate paper target unless narrowed.

### Better narrowed versions
- periodic or switched defense schedules
- actuator saturation under destabilizing perturbations
- packet-drop attack channels with stochastic losses

### Attempted solution path
- start with periodic switching between a small set of masks
- study averaged or lifted dynamics
- compare with simultaneous mixed strategies

### Code note
Only pursue after the LTI finite-action baseline is finished.

---

## 2.14 The naive full-FMP “hide or rotate the dynamics” question
**Priority:** 14  
**Assessment:** weak as originally posed  
**Main source questions covered:** the first paragraph of 4.1.3

### Why it is weak
The future-work chapter already diagnoses the issue correctly: if you know where the attacker will strike and can arbitrarily redesign the realization, then “hiding” vulnerable dynamics from those coordinates can trivialize the problem. That makes the result overly dependent on an unrealistic fixed exposure model.

### Recommendation
Do **not** make this a main research thrust. Keep it only as:
- a motivating failure mode,
- a counterexample to a naive formulation,
- a reason to adopt the random-exposure expected-vulnerability reformulation instead.

### Code note
Only include a tiny demonstration showing how the naive problem can produce misleading “security by rotation” artifacts.

---

## 2.15 Nonlinear ISS / small-gain masked perturbations
**Priority:** 15  
**Assessment:** good, but adjacent and longer-horizon  
**Source:** our discussion, not only the PDF

### Why include it
It is not one of the explicit Chapter 4 questions, but it is a very natural extension of the same program and likely a stronger nonlinear story than “attacks designed from linearization destabilize the nonlinear plant.”

### Core idea
Model:
- plant as an IOS/ISS system from `(u,w)` to `(x,r)`
- attacker as an IOS/ISS system from `r` to `w`
- define a safe set of attacks through a nonlinear small-gain condition

### Better formulation than the original sketch
The key condition should not merely be “the attacker’s `w` ultimate bound lies inside a predeclared box,” but a **self-consistent loop gain condition** such as:

$$
\sigma_w \circ \rho_\Delta(s) < s
$$

on a relevant domain, where:
- `σ_w` is the plant-side IOS gain from `w` to `r`
- `ρ_Δ` is the attacker-side IOS gain from `r` to `w`

### Assessment
Very good theoretically, but not the fastest next project if the masked-LTI line still has many publishable open questions.

### Recommended background
- ISS / IOS small-gain theorems
- nonlinear interconnected-system stability
- safety verification via small gain
- Khalil; Jiang, Teel, Praly; Dashkovskiy-Rüffer-Wirth

### Code note
Only a toy simulation or theorem note is appropriate at first.

---

# 3. Recommended research order by quarter

## Phase A: immediate, most likely to produce usable results
1. Simultaneous destabilization numerics and support-structure lemmas
2. Reduced mixed-defense action-set algorithms
3. Stackelberg variant on finite reduced games
4. Signal-importance heuristics
5. Budget-monotonicity/counterexample search

## Phase B: second-wave theory and stochastic reformulations
6. Expected-vulnerability full-FMP reformulation
7. Better defender objectives and Pareto frontiers
8. Forcing attacker mixing
9. Parametric sensitivity / policy-region maps

## Phase C: harder theory
10. Interpolation certificates and impossibility results
11. Partial/probabilistic masks as uncertainty sets
12. Nonlinear ISS extension
13. Time-varying and non-LTI channels

---

# 4. Recommended paper/dissertation packaging

## Paper 1: simultaneous destabilization and attack dominance
- support-overlap lemmas
- pairwise simultaneous destabilization numerics
- dominant attack constructions and counterexamples

## Paper 2: reduced mixed defense via success-set geometry
- defense-subset selection
- approximation algorithms
- reduced zero-sum / Stackelberg evaluation

## Paper 3: full FMP under exposure uncertainty
- why the naive formulation fails
- stochastic exposure reformulation
- toy realization-design experiments

## Paper 4: nonlinear small-gain extension
- safe-set definition via IOS/ISS gains
- theorem and toy simulations

---

# 5. Suggested bibliography to prioritize

## Core robust-control / vulnerability background
- Zhou, Doyle, Glover — *Robust and Optimal Control*
- Basar & Bernhard — *H-infinity Optimal Control and Related Minimax Design Problems*
- Rai, Ward, Roy, Warnick (2012) — vulnerable links / secure architectures
- Cox, Roy, Warnick (2014) — science of system security
- Chetty et al. (2014) — distributed/coordinated destabilization attacks
- Chetty & Warnick (2020) — structure in networks of dynamic systems

## Simultaneous stabilization / interpolation
- Blondel — *Simultaneous stabilization of linear systems*
- Doyle, Francis, Tannenbaum — *Feedback Control Theory*
- Cui & Lindquist (2023, 2024) — simultaneous stabilization by analytic interpolation
- literature on strong stabilization and Nevanlinna–Pick interpolation

## Security games / approximation
- Kiekintveld et al. on randomized resource allocation for security games
- Sinha, Fang, An, Kiekintveld, Tambe — Stackelberg security games survey
- Balcan et al. — commitment without regret
- set cover / max coverage / submodular maximization classics

## Full-FMP / realization / DSF
- DSF identification / structure literature already cited in your thesis
- structured realization papers
- stochastic / risk-sensitive design references

## Nonlinear extension
- Khalil — *Nonlinear Systems*
- Jiang, Teel, Praly on nonlinear small gain
- Dashkovskiy, Rüffer, Wirth on ISS small-gain for large-scale systems
- more recent small-gain/safety verification papers

---

# 6. One consolidated Codex prompt for the next coding wave

```text
Extend the masked-perturbation Python package with the next-wave research modules.

This wave should focus on:
1. simultaneous destabilization across multiple defense masks,
2. reduced mixed-defense action-set design from success sets,
3. Stackelberg / defender-commitment games,
4. signal-importance heuristics,
5. budget-monotonicity counterexample search,
6. optional toy expected-vulnerability studies for the full FMP.

General requirements:
- keep all functionality modular and documented
- support small exact examples first
- never claim approximation guarantees unless they are actually justified
- produce useful plots, markdown reports, and unit tests
- preserve the transfer-function / masked-map formulation as the public abstraction

Modules to add:
- simultaneous.py
- defense_reduction.py
- stackelberg.py
- signal_importance.py
- counterexamples.py
- stochastic_exposure.py
- reports_future_work.py

Required features:

A. Simultaneous destabilization
- input: base model map M, finite mask set {∇_j}, finite admissible attack family or attack generators
- compute candidate minimal attacks per mask
- test pairwise and multi-mask simultaneous success
- test additive combined attacks
- build support-overlap and cross-success matrices
- output coverage tables and figures

B. Reduced defense action-set design
- compute success sets S(∇)
- implement greedy coverage, greedy overlap-based selection, beam search, local search
- solve reduced zero-sum game on the chosen defense subset
- compare against full game on small instances
- plot value vs subset size z

C. Stackelberg game
- implement defender commitment optimization with strong Stackelberg and adversarial tie-breaking variants
- compare Stackelberg value to simultaneous minimax value
- integrate with dominated-strategy elimination

D. Signal importance
- marginal vulnerability-drop scoring
- success-set marginal coverage scoring
- optional pairwise / Monte Carlo Shapley approximations
- greedy defense construction from ranked signals

E. Budget monotonicity / counterexamples
- sweep defender budget
- detect instances where more budget does not help or hurts
- save minimal counterexamples and produce a markdown report

F. Stochastic exposure / toy full-FMP study
- sample exposure subsets from a user-specified distribution
- compute expected and worst-case vulnerability across exposures
- compare designs under fixed known exposure vs random exposure

Outputs:
- notebooks:
  - simultaneous_destabilization_demo.ipynb
  - defense_subset_design_demo.ipynb
  - stackelberg_demo.ipynb
  - signal_importance_demo.ipynb
  - budget_counterexamples_demo.ipynb
  - stochastic_exposure_demo.ipynb
- markdown report summarizing all experiments
- unit tests for all exact small examples

Implementation notes:
- use scipy / numpy / python-control / matplotlib / pandas
- keep examples small and interpretable
- prefer exact exhaustive validation whenever state spaces are tiny
- tag each result as exact, heuristic, or exploratory
```

---

# 7. Short recommendations

## Best next theorem-driven project
**Simultaneous destabilization with support-structure lemmas and pairwise constructions.**

## Best next algorithmic project
**Mixed-defense action-set reduction via success-set geometry and approximation heuristics.**

## Best next practical/game-theoretic project
**Stackelberg masked defense with reduced action sets.**

## Best “do not overinvest in this exact form”
**Naive full FMP by hiding/rotating dynamics away from fixed known exposed signals.**

## Best long-horizon extension
**Nonlinear ISS / small-gain safe-set extension.**
