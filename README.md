# Flu Virus Spread in a City

A GIS-based agent simulation of influenza transmission in a city, implemented
in GAML for GAMA Platform 2025.6.

The model combines population movement, household and workplace contacts, the
SIR epidemic model, local testing and isolation, vaccination, virus mutation,
variant-specific immunity, and batch experimentation.

The main research question is:

> How can local isolation reduce the diffusion of an epidemic?

## Table of contents

- [Project structure](#project-structure)
- [Requirements implemented](#requirements-implemented)
- [How to run the model](#how-to-run-the-model)
- [Model architecture](#model-architecture)
- [Global configuration](#global-configuration)
- [Initialization](#initialization)
- [Daily public-health process](#daily-public-health-process)
- [Agent species](#agent-species)
- [Movement and animation](#movement-and-animation)
- [Disease transmission](#disease-transmission)
- [Recovery and immunity](#recovery-and-immunity)
- [Vaccination](#vaccination)
- [Virus variants](#virus-variants)
- [Simulation timeline](#simulation-timeline)
- [Indicators and outputs](#indicators-and-outputs)
- [Experiments](#experiments)
- [How to reproduce model versions](#how-to-reproduce-model-versions)
- [Assumptions and interpretation](#assumptions-and-interpretation)
- [Extending the model](#extending-the-model)
- [Validation](#validation)

## Project structure

```text
FluVirus/
|- models/
|  `- FluVirus.gaml
|- includes/
|  |- people.obj
|  |- building_environment.shp
|  |- building_environment.dbf
|  |- building_environment.shx
|  |- road_environment.shp
|  |- road_environment.dbf
|  |- road_environment.shx
|  `- other supporting GIS files
|- .project
`- README.md
```

The `.shp`, `.dbf`, `.shx`, and related files are parts of the same GIS
datasets. Keep them together in `includes/`; moving only the `.shp` file will
usually prevent the GIS layer from loading correctly.

## Requirements implemented

The model implements the project specification as follows:

| Requirement | Implementation |
|---|---|
| GIS city | Buildings and roads are loaded from shapefiles. |
| Population | Every residential building receives one family of 3-6 people. |
| Children | Each family has 0-2 children. |
| School | The largest building is designated as the school. |
| Workplaces | 25% of non-school buildings are workplaces by default; each adult receives one. |
| Daily mobility | People leave at 08:00 and return home at 18:00. |
| Road movement | People use GAMA's `moving` skill and the GIS road graph. |
| 3D city | GIS buildings are extruded by role and rendered with lighting in a native 3D display. |
| SIR states | People are susceptible, infected, or recovered. |
| Infection | Transmission probability is 33% per evaluated contact, with at most 10 contacts per event by default. |
| Recovery | Infected people recover after 10 simulated days. |
| Testing | Authorities test 1% of the population each day. |
| Isolation | A detected infected person stays home for 12 days. |
| Vaccination | Authorities vaccinate 0.05% of the population per day. |
| Vaccine effect | Matched vaccination divides transmission risk by 3. |
| Variants | Infections have a daily 0.1% mutation probability. |
| Variant transmissibility | A new variant uses parent probability times a random value from 0.5 to 1.5. |
| Variant vaccine escape | Protection against an unmatched variant is two times weaker. |
| Reinfection | A recovered person can be infected by a different variant. |
| Batch analysis | Starting vaccination coverage is tested at 10%, 50%, and 90%. |

## How to run the model

### Interactive simulation

1. Open GAMA Platform 2025.6.
2. Import or open the `FluVirus` project directory.
3. Open `models/FluVirus.gaml`.
4. Select the `flu_city` experiment.
5. Click Run.
6. Use the experiment parameters to enable or disable interventions.

The experiment opens two main views:

- `3D Flu City` shows a landscaped city under a blue sky, moving clouds, a sun,
  detailed homes, workplaces, the school, roads, and animated people.
- `Epidemic dashboard` shows SIR curves and public-health intervention curves.

The 3D city uses native GAMA rendering rather than an external graphics
library. Drag to rotate the scene, scroll to zoom, and use the display toolbar
to reset the camera when needed.

With **Show detailed 3D buildings** enabled, homes receive pitched orange
roofs, glowing windows, doors, and chimneys; workplaces receive glass window
bands, entrance canopies, rooftop service blocks, and antennas; and the school
receives a columned entrance, flagpole, and blue-roofed central clock tower.
The **Show sky and moving clouds** option adds a blue atmosphere, stylized sun,
and eleven softly shaded cloud groups that drift slowly across the city. With
**Use realistic 3D people** enabled, each person is rendered with the human mesh
in `includes/people.obj`.
Adults and children use different scales, agents rotate toward their movement
heading, and clothing receives an individual color. These decorations are
visual only and do not change the GIS footprint, destination, or epidemic
logic.

### Native GAMA examples used for the 3D design

The implementation follows reusable patterns from the open model library
bundled with GAMA 2025.6.4:

- **Building Elevation** — extrudes building shapefile polygons with the
  `depth:` facet and uses a lit 3D display;
- **Luneray Flu, Model 6** — combines 3D GIS buildings, road geometry, and
  animated OBJ human models in a flu simulation;
- **Incremental Model 7** — renders health-colored people as spheres above
  extruded buildings;
- **Road Traffic, Movement of People** — constrains `goto` movement with
  `on: road_network`.

These examples are available inside GAMA under **Library models**. The project
packages its human OBJ mesh in `includes/people.obj` and otherwise uses native
GAML shapes and lighting, so no extra plugin or download is required.

### Automated test

Run the `automated_smoke_test` batch experiment to verify that the local GAMA
installation can:

- load the GIS files;
- build the road network;
- create all families;
- perform repeated home, work, and school trips;
- execute infection, testing, isolation, vaccination, and mutation logic;
- complete the 10-day recovery lifecycle.

This test runs for 11 simulated days without opening visualization displays.

## Model architecture

The model contains one global controller and three agent species:

```text
global model
|- loads GIS data
|- creates the road graph
|- creates families
|- manages testing, vaccination, and variants
|- calculates global indicators
|
|- road
|  `- represents a segment of the movement network
|
|- building
|  `- represents homes, workplaces, and the school
|
`- person
   |- stores health and vaccination state
   |- travels between buildings
   |- transmits infection
   |- recovers and stores variant immunity
   `- follows isolation rules
```

The source begins with:

```gaml
model FluVirus
```

This declares the GAML model name. The remainder of the file defines its
global state, species, actions, reflexes, aspects, and experiments.

## Global configuration

### Simulation clock

```gaml
float step <- 10#mn;
int cycles_per_hour <- 6;
int cycles_per_day <- 144;
```

One simulation cycle represents ten minutes. Therefore:

- 6 cycles = 1 hour;
- 144 cycles = 1 day;
- 1,440 cycles = 10 days.

Ten-minute cycles allow road movement to remain visible over several frames.
A one-hour cycle would make many trips appear instantaneous on this compact
GIS map.

The displayed date and time are calculated from the cycle:

```gaml
int simulation_day -> int(cycle / cycles_per_day);
int hour_of_day -> int((cycle mod cycles_per_day) / cycles_per_hour);
int minute_of_hour -> int((cycle mod cycles_per_hour) * 10);
```

The `->` operator means the value is recalculated whenever it is requested.

### GIS inputs

```gaml
file road_file <- file("../includes/road_environment.shp");
file building_file <- file("../includes/building_environment.shp");
geometry shape <- envelope(road_file) + envelope(building_file);
graph<geometry, geometry> road_network;
building school;
```

- `road_file` contains road line geometries.
- `building_file` contains building polygons.
- `shape` defines the simulation boundary from both GIS envelopes.
- `road_network` stores a graph created from road geometries.
- `school` stores a reference to the selected school building.

### Epidemiological parameters

| Variable | Default | Meaning |
|---|---:|---|
| `base_infection_probability` | `0.33` | Infection probability per evaluated contact for the original variant. |
| `infectious_period_days` | `10` | Days from infection until recovery. |
| `daily_testing_rate` | `0.01` | Proportion of the population randomly tested each day. |
| `isolation_period_days` | `12` | Duration of home isolation after detection. |
| `initially_infected` | `8` | Number of infected agents at initialization. |

### Extension controls

| Variable | Default | Meaning |
|---|---:|---|
| `enable_isolation` | `true` | Enables testing and home isolation. |
| `enable_vaccination` | `true` | Enables initial and daily vaccination. |
| `enable_variants` | `true` | Enables mutation and vaccine retargeting. |
| `initial_vaccination_coverage` | `0.10` | Starting vaccinated fraction. |
| `daily_vaccination_rate` | `0.0005` | Daily vaccination capacity, equal to 0.05%. |
| `daily_mutation_probability` | `0.001` | Mutation probability per infected person per day. |

### Global indicators

The following dynamic values count the current population states:

```gaml
int population_size -> length(person);
int susceptible_count -> length(person where (each.health_state = 0));
int infected_count -> length(person where (each.health_state = 1));
int recovered_count -> length(person where (each.health_state = 2));
int isolated_count -> length(person where each.isolated);
int vaccinated_count -> length(person where each.vaccinated);
int commuting_count -> length(person where each.travelling);
```

Additional indicators include:

- `vaccinated_percent`: current vaccinated share of the population;
- `cumulative_infections`: total infection events, including reinfections;
- `attack_rate`: cumulative infection events divided by population size;
- `peak_infected`: largest simultaneous infected count observed;
- `vaccine_target_variant`: variant currently targeted by vaccination;
- `epidemic_finished`: true after day zero when no infected person remains.

Because `cumulative_infections` includes reinfections, `attack_rate` is an
event-based measure and can exceed 100% if many people are infected more than
once by different variants.

## Initialization

The global `init` block constructs the complete simulation.

### 1. Load GIS agents

```gaml
create road from: road_file;
create building from: building_file;
road_network <- as_edge_graph(road);
```

Each road geometry becomes a `road` agent and each building polygon becomes a
`building` agent. `as_edge_graph(road)` converts the road agents into a graph
that the moving skill can use for pathfinding.

### 2. Select the school

```gaml
school <- building with_max_of (each.shape.area);
ask school {
    is_school <- true;
}
```

The building with the largest polygon area is marked as the school.

### 3. Assign homes and workplaces, then create families

After reserving the school, the model selects the configured share of the
remaining buildings as workplaces. All remaining buildings become residential
homes. For every residential building, the model randomly chooses:

- a family size from 3 to 6;
- a child count from 0 to 2, while ensuring at least one adult.

Children use the school as their `day_place`. Adults receive a random building
from the exclusive workplace list.
Every person starts at a random point inside the home polygon.

The total population is not entered directly. It is produced from the number
of GIS buildings and their randomly generated family sizes.

### 4. Apply initial vaccination

If vaccination is enabled, the model randomly selects the requested starting
percentage of people and marks them as vaccinated against variant 0.

### 5. Seed infections

The model randomly selects up to `initially_infected` people and calls:

```gaml
do become_infected(0, base_infection_probability);
```

Variant 0 is the original virus strain.

## Daily public-health process

The global `daily_public_health` reflex executes once every day at 06:00.

```gaml
reflex daily_public_health when: (cycle > 0) and
    ((cycle mod cycles_per_day) = (6 * cycles_per_hour))
```

It performs three processes in order.

### Testing and isolation

The number tested is:

```text
max(1, integer(population size x daily testing rate))
```

The selected people are random, so testing represents population sampling
rather than contact tracing or symptom-based testing.

If a tested person is infected and not already isolated:

1. `isolated` becomes true;
2. `isolation_end_day` is set to the current day plus 12;
3. `return_home` is called immediately.

An isolated person cannot begin the morning trip to work or school.

### Daily vaccination

At the specified rate, a population near 1,000 produces about half a dose per
day. Rounding that value every day would incorrectly create either zero or one
full dose. The model therefore uses `vaccination_credit` as an accumulator.

Example for 0.5 doses per day:

```text
Day 1: credit = 0.5, doses = 0, remaining credit = 0.5
Day 2: credit = 1.0, doses = 1, remaining credit = 0.0
Day 3: credit = 0.5, doses = 0, remaining credit = 0.5
```

This preserves the requested average rate over time. Doses are given only to
unvaccinated people and target `vaccine_target_variant`.

### Mutation and vaccine surveillance

Each infected person receives an independent daily mutation trial. When a
mutation occurs:

1. `next_variant_id` is incremented;
2. the infection receives the new variant ID;
3. its transmission probability becomes the parent probability multiplied by
   a random value between 0.5 and 1.5;
4. the resulting probability is limited to the range 0.01-0.95.

The model then counts active cases of each circulating variant. The most
prevalent variant becomes the new vaccine target.

## Agent species

### `road`

Road agents are passive GIS objects. Their aspect draws each geometry in gray:

```gaml
draw shape color: rgb(108, 117, 125) width: 2.2#m;
```

Their main functional role is to form `road_network` for person pathfinding.

### `building`

Each building stores:

```gaml
bool is_school <- false;
bool is_home <- false;
bool is_workplace <- false;
float display_height <- 8#m;
```

Building roles are exclusive. The largest GIS polygon becomes the school,
25% of the remaining buildings are workplaces by default, and all other
buildings are family homes. Homes are cyan, workplaces are pink, and the
school is gold. In the 3D aspect, every GIS polygon is extruded by
`display_height`; school and workplace roofs also receive labels.

### `person`

The `person` species includes the GAMA moving skill:

```gaml
species person skills: [moving]
```

The skill provides movement variables and actions such as `speed`,
`current_path`, and `goto`.

#### Health fields

| Field | Purpose |
|---|---|
| `health_state` | `0` susceptible, `1` infected, `2` recovered. |
| `infection_start_day` | Day the current infection began. |
| `active_variant` | Current infecting variant, or `-1` if not infected. |
| `infection_probability` | Transmission probability of the active variant. |
| `immune_variants` | Variants from which the person has recovered. |

#### Demographic and location fields

| Field | Purpose |
|---|---|
| `is_child` | Distinguishes children from adults. |
| `home` | Family building. |
| `day_place` | School for a child or workplace for an adult. |
| `current_building` | Building currently occupied; `nil` while travelling. |
| `target_building` | Destination building for the active trip. |
| `travel_target` | Exact point inside the destination building. |
| `travelling` | Whether movement is currently active. |

#### Intervention fields

| Field | Purpose |
|---|---|
| `isolated` | Whether the person must remain at home. |
| `isolation_end_day` | Day isolation finishes. |
| `vaccinated` | Whether any vaccine has been received. |
| `vaccine_variant` | Variant targeted by the person's vaccine. |

## Movement and animation

Movement follows the road-network pattern used in the reference evacuation
model.

### Starting a trip

`start_trip` records the target building, chooses a precise point inside it,
clears any previous path, removes the person from the current building, and
sets `travelling` to true.

```gaml
action start_trip (building destination_building) {
    target_building <- destination_building;
    travel_target <- any_location_in(destination_building.shape);
    current_path <- nil;
    current_building <- nil;
    travelling <- true;
}
```

### Following roads

While travelling, the agent executes:

```gaml
do goto target: travel_target on: road_network recompute_path: false;
```

Important parts:

- `target` is an exact point in the destination building;
- `on: road_network` constrains the route to the GIS road graph;
- `recompute_path: false` reuses the route instead of calculating a new
  shortest path on every animation frame.

When the agent reaches `travel_target`, the model registers arrival, clears the
destination and cached path, and stops travel.

### Daily travel schedule

| Time | Behavior |
|---|---|
| 08:00 | Non-isolated adults travel to work and children travel to school. |
| During trip | People move along the road graph for multiple display frames. |
| 18:00 | People begin the return trip home. |
| Any time after isolation detection | The detected person begins returning home. |

Children move at 0.8 km/h and adults at a random 0.9-1.3 km/h. These speeds are
scaled for the compact 500 x 450 GIS environment so that movement remains
visible and does not look like teleportation.

## Disease transmission

Contacts are evaluated twice per day:

- 12:00, representing daytime workplace and school contact;
- 22:00, representing evening household contact.

An infected person selects potential contacts who:

- are not the infected person;
- are not currently travelling;
- occupy the same building;
- are not already infected;
- are susceptible, or recovered without immunity to the active variant.

Each infectious agent meets at most `max_contacts_per_event` eligible people
at an event. This prevents a large school from behaving as if every child had
close contact with every other child every day. For each sampled contact, the
model starts with the infecting person's variant transmission probability.

### Vaccine-adjusted probability

For an unvaccinated contact:

```text
effective probability = variant transmission probability
```

For a contact vaccinated against the same variant:

```text
effective probability = variant transmission probability / 3
```

For a vaccinated contact facing a different variant:

```text
effective probability = variant transmission probability / 1.5
```

Dividing by 1.5 instead of 3 makes the unmatched vaccine two times less
effective than the matched vaccine.

If the Bernoulli trial succeeds, `become_infected` records the variant,
probability, start day, and increments `cumulative_infections`.

## Recovery and immunity

The recovery reflex checks:

```gaml
simulation_day - infection_start_day >= infectious_period_days
```

After the default 10 days:

1. health state changes to recovered;
2. the recovered variant is added to `immune_variants`;
3. `active_variant` is cleared.

Recovery protects against the same variant but not against a variant that is
not present in `immune_variants`. This implements variant-specific reinfection.

Isolation and infection use separate clocks. A person can recover after 10
days but remain isolated until the 12-day isolation period ends.

## Vaccination

Vaccination does not create a separate SIR health state. A vaccinated person
can still be susceptible, infected, or recovered.

The vaccine changes the probability of infection; it does not guarantee
complete immunity. The white outline on the map shows vaccinated people.

When authorities update `vaccine_target_variant`, only future vaccinations use
the new target. Existing vaccinated people keep the variant ID of the vaccine
they previously received.

## Virus variants

Variant 0 is the original strain. New variants receive IDs 1, 2, 3, and so on.

Each infection stores its own:

- variant ID;
- transmission probability.

This allows descendants of different infection chains to have different
transmissibility. The mutation bounds prevent impossible negative probabilities
and cap extremely transmissible variants at 95% per evaluated contact.

## Simulation timeline

A normal simulated day follows this schedule:

```text
06:00  Test a random 1% sample
       Isolate detected infected people
       Allocate accumulated vaccination doses
       Mutate active infections where applicable
       Select the most prevalent variant as vaccine target

08:00  Adults depart for workplaces
       Children depart for school

08:00-  Agents animate along the road network until arrival

12:00  Evaluate workplace and school transmission

18:00  Everyone not already home starts the return trip

22:00  Evaluate household transmission

All day Recover infections that reached the infectious-period limit
        Release people whose isolation period has ended
        Update peak infection statistics
```

## Indicators and outputs

### Monitors

The interactive experiment shows:

| Monitor | Meaning |
|---|---|
| Day / time | Current simulated day and clock time. |
| Susceptible | Number currently in S state. |
| Infected | Number currently in I state. |
| Recovered | Number currently in R state. |
| Isolated | Number under home isolation. |
| Commuting on roads | Number currently following a road route. |
| Homes / workplaces / schools | Number of GIS buildings assigned to each exclusive role. |
| Vaccinated (%) | Percentage ever vaccinated. |
| Peak infected | Maximum simultaneous infections observed. |
| Cumulative attack rate (%) | Infection events divided by population size. |
| Current vaccine target | Variant targeted by new doses. |

### Map legend

| Appearance | Meaning |
|---|---|
| Green beacon above person | Susceptible. |
| Red beacon above person | Infected. |
| Blue beacon above person | Recovered. |
| Purple outline | Isolated. |
| White outline | Vaccinated and not isolated. |
| Smaller human model | Child. |
| Larger human model | Adult. |
| Gold building with blue clock tower | School. |
| Cyan building with orange pitched roof | Family home. |
| Pink tower with rooftop block and antenna | Workplace. |
| Dark gray line | Road. |
| Blue background, yellow sun, white clouds | Procedural sky and atmosphere. |

Isolation has visual priority over vaccination: an isolated vaccinated person
uses the purple isolation outline.

### Charts

`SIR population through time` plots susceptible, infected, and recovered
counts. `Public-health response` plots isolated and vaccinated counts.

Charts refresh every six cycles, which equals once per simulated hour.

## Experiments

### `flu_city`

This is the main GUI experiment. It exposes parameters in five groups:

- Epidemic;
- City roles;
- Public health;
- Vaccination;
- Variants.

Changing parameters and restarting the experiment makes it possible to compare
policy scenarios without editing the source code.

### `vaccination_coverage_analysis`

This batch experiment evaluates starting vaccination coverage values:

```text
10%, 50%, 90%
```

Each coverage is repeated 10 times with controlled seeds. A run stops when the
epidemic finishes or reaches day 120. Testing, vaccination, and variants are
enabled for all scenarios.

In the GAMA batch interface, compare the final and peak indicators among the
three coverage settings. The model does not currently write a CSV file; add a
`save` statement if external statistical analysis is required.

### `automated_smoke_test`

This non-GUI batch experiment runs one 11-day simulation. It is intended for
compilation and runtime validation, not policy analysis.

## How to reproduce model versions

The incremental versions from the presentation can be reproduced from the
`flu_city` parameter panel.

### Version 1 - Baseline SIR

Set:

```text
Enable testing and isolation = false
Enable vaccination = false
Enable variants = false
```

This leaves GIS movement, building contacts, infection, and recovery.

### Version 2 - Local isolation

Set:

```text
Enable testing and isolation = true
Enable vaccination = false
Enable variants = false
```

Compare peak infections and attack rate with Version 1.

### Version 3 - GIS structure

GIS buildings, road movement, and spatial destinations are always active in
this implementation.

### Version 4 - Families and school

Family creation, children, household assignment, and school travel are always
active in this implementation.

### Extensions

Enable vaccination and variants individually or together. Use the batch
experiment to compare initial vaccination coverage.

## Assumptions and interpretation

The implementation makes several explicit modeling assumptions:

1. Each infected agent meets a configurable random sample of at most 10
   eligible people at each scheduled contact event by default.
2. Transmission is evaluated at noon and 22:00, not continuously every cycle.
3. Testing is random and has perfect sensitivity and specificity.
4. A positive test causes immediate compliance with isolation.
5. Isolation lasts 12 days even if recovery happens earlier.
6. Recovery occurs after a fixed 10-day infectious period.
7. Recovered immunity is variant-specific.
8. Vaccination reduces infection probability but does not affect recovery time.
9. Each person receives at most one vaccine in the current implementation.
10. The vaccine target can change, but previously vaccinated people are not
    automatically revaccinated.
11. Mutation occurs within an infected person and replaces that person's active
    variant.
12. The largest GIS building is assumed to be a school.
13. Roads constrain movement, but congestion and vehicle interactions are not
    modeled.
14. Movement speeds are visualization-scaled to this GIS dataset.

These assumptions should be stated when interpreting results. The simulation
is a research and teaching model, not a clinical forecasting system.

## Extending the model

Possible next improvements include:

- export daily indicators and batch results to CSV;
- add age-dependent infection and recovery probabilities;
- distinguish household, school, workplace, and road transmission rates;
- add test sensitivity, specificity, and delayed results;
- allow booster vaccination and revaccination against a new target variant;
- add symptom severity and asymptomatic infections;
- create workplace capacities instead of selecting workplaces uniformly;
- model road congestion or public transport exposure;
- add deaths, hospitalization, or healthcare capacity;
- run formal sensitivity analysis over testing rate, isolation duration,
  transmission probability, and vaccination coverage.

## Validation

The final model was checked with GAMA Platform 2025.6.4.

Validation completed:

- GAML compilation succeeded;
- GIS road and building files loaded successfully;
- the road graph was created successfully;
- families, children, workplaces, and the school were created;
- the animated morning and evening road trips executed without runtime errors;
- cached road paths worked with `recompute_path: false`;
- daily testing, isolation, vaccination, and mutation executed;
- infection and variant-specific reinfection logic executed;
- the 10-day recovery lifecycle completed in the 11-day smoke test.

The main implementation is located at `models/FluVirus.gaml`.
