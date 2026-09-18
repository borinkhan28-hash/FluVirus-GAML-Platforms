/**
 * Flu Virus Spread in a City
 *
 * GIS-based SIR model implementing the USTH Project 5 specification:
 * commuting, 33% contact transmission, 10-day recovery, testing and
 * 12-day local isolation, families, a school, vaccination, variants,
 * and a vaccination-coverage batch experiment.
 */
model FluVirus

global {
	// Faster clock: one cycle is 30 simulated minutes instead of 10, so a day
	// takes 48 cycles instead of 144 (3x fewer steps for the same number of days).
	// Edit minutes_per_step only; the derived values follow automatically.
	int minutes_per_step <- 30;
	float step <- minutes_per_step * 1#mn;
	int cycles_per_hour <- 60 / minutes_per_step;
	int cycles_per_day <- 24 * cycles_per_hour;

	// GIS inputs bundled with this project.
	file road_file <- file("../includes/road_environment.shp");
	file building_file <- file("../includes/building_environment.shp");
	geometry shape <- envelope(road_file) + envelope(building_file);
	graph<geometry, geometry> road_network;
	building school;
	list<building> home_buildings <- [];
	list<building> workplace_buildings <- [];
	float workplace_building_fraction <- 0.25;
	bool show_cartoon_details <- true;
	bool use_realistic_people_3d <- true;
	bool show_sky_and_clouds <- true;

	// Core epidemiological assumptions from the assignment.	
	float base_infection_probability <- 0.33;

	int infectious_period_days <- 10;
	float daily_testing_rate <- 0.01;
	int isolation_period_days <- 12;
	int initially_infected <- 8;
	int max_contacts_per_event <- 10;
	int departure_hour <- 8;
	int return_hour <- 18;
	int work_contact_hour <- 12;
	int home_contact_hour <- 22;

	// Optional extensions. They can be switched off for a V1/V2 baseline.
	bool enable_isolation <- true;
	bool enable_vaccination <- true;
	bool enable_variants <- true;
	// City-structure extensions. Disabling families spawns one adult per home;
	// disabling the school also removes children, since they have no day_place.
	bool enable_families <- true;
	bool enable_school <- true;
	// The GIS city is two blocks joined by a single 125 m road across the
	// central gap. Disabling this connector splits the road graph in two and
	// keeps the communities epidemiologically separate.
	bool enable_block_connector <- true;
	float initial_vaccination_coverage <- 0.10;
	float daily_vaccination_rate <- 0.0005;
	float daily_mutation_probability <- 0.001;
	// When enabled, the interactive experiment stops as soon as the last
	// infection clears (epidemic_finished).
	bool stop_when_epidemic_finished <- false;

	// Runtime state and indicators.
	int simulation_day -> int(cycle / cycles_per_day);
	float simulation_time_days -> float(cycle) / cycles_per_day;
	int hour_of_day -> int((cycle mod cycles_per_day) / cycles_per_hour);
	int minute_of_hour -> int((cycle mod cycles_per_hour) * minutes_per_step);
	int population_size -> length(person);
	int susceptible_count -> length(person where (each.health_state = 0));
	int infected_count -> length(person where (each.health_state = 1));
	int recovered_count -> length(person where (each.health_state = 2));
	int isolated_count -> length(person where each.isolated);
	int vaccinated_count -> length(person where each.vaccinated);
	int commuting_count -> length(person where each.travelling);
	float vaccinated_percent -> population_size = 0 ? 0.0 : 100.0 * vaccinated_count / population_size;
	float attack_rate -> population_size = 0 ? 0.0 : 100.0 * cumulative_infections / population_size;
	int cumulative_infections <- 0;
	int peak_infected <- 0;
	int next_variant_id <- 0;
	int vaccine_target_variant <- 0;
	float vaccination_credit <- 0.0;
	bool epidemic_finished -> (simulation_day > 0) and (infected_count = 0);

	init {
		create road from: road_file;
		create building from: building_file;

		// The connector is the only road that crosses the central gap, so it is
		// the one segment taller than the 25 m city grid. Flag it, then leave it
		// out of the routing graph when the connector is disabled.
		ask road {
			is_connector <- (shape.height > 100#m) and (shape.width < 100#m);
		}
		road_network <- enable_block_connector
			? as_edge_graph(road)
			: as_edge_graph(road where !each.is_connector);
		create sky_cloud number: 11 {
			location <- any_location_in(shape);
			cloud_altitude <- rnd(48#m, 76#m);
			cloud_size <- rnd(10#m, 18#m);
			drift_speed <- rnd(0.025#m, 0.065#m);
		}

		// Building-role logic follows the separated-role reference model:
		// one school, a configurable workplace group, and family homes.
		if enable_school {
			school <- building with_max_of (each.shape.area);
			ask school {
				is_school <- true;
				display_height <- 20#m;
			}
		}
		list<building> available_buildings <- (school = nil) ? building
			: (building where (each != school));
		int workplace_count <- int(length(available_buildings) * workplace_building_fraction);
		workplace_count <- max([1, workplace_count]);
		workplace_count <- min([workplace_count, length(available_buildings) - 1]);
		workplace_buildings <- workplace_count among available_buildings;
		home_buildings <- available_buildings - workplace_buildings;
		ask workplace_buildings {
			is_workplace <- true;
			display_height <- rnd(13#m, 22#m);
		}
		ask home_buildings {
			is_home <- true;
			display_height <- rnd(6#m, 12#m);
		}

		// Split the city into two communities at the central gap (y 150-275),
		// matching the two road blocks joined only by the connector road.
		ask building {
			community <- (location.y > 212#m) ? 1 : 0;
		}

		// Every residential building receives a household. With families on
		// this is 3-6 people (0-2 children); adults work and children attend
		// the school. With families off, one independent adult lives per home.
		ask home_buildings {
			int family_size <- rnd(3, 6);
			int child_count <- (enable_families and enable_school)
				? rnd(0, min([2, family_size - 1])) : 0;

			if enable_families {
				create person number: child_count {
					home <- myself;
					is_child <- true;
					day_place <- school;
					current_building <- home;
					location <- any_location_in(home.shape);
				}

				create person number: family_size - child_count {
					home <- myself;
					is_child <- false;
					day_place <- one_of(workplace_buildings);
					current_building <- home;
					location <- any_location_in(home.shape);
				}
			} else {
				create person number: 1 {
					home <- myself;
					is_child <- false;
					day_place <- one_of(workplace_buildings);
					current_building <- home;
					location <- any_location_in(home.shape);
				}
			}
		}

		// Each person belongs to the community of their home building.
		ask person {
			community <- home.community;
		}

		// Apply starting coverage first, then seed infections among everybody.
		if enable_vaccination and (initial_vaccination_coverage > 0.0) {
			int initial_doses <- min([population_size, int(population_size * initial_vaccination_coverage)]);
			ask initial_doses among person {
				vaccinated <- true;
				vaccine_variant <- 0;
			}
		}

		int seed_cases <- min([initially_infected, population_size]);
		ask seed_cases among person {
			do become_infected(0, base_infection_probability);
		}
	}

	// Testing, isolation, vaccination, variant surveillance and mutation occur daily.
	reflex daily_public_health when: (cycle > 0) and
		((cycle mod cycles_per_day) = (6 * cycles_per_hour)) {
		if enable_isolation {
			int tests_today <- min([population_size, max([1, int(population_size * daily_testing_rate)])]);
			ask tests_today among person {
				if (health_state = 1) and !isolated {
					isolated <- true;
					isolation_end_day <- simulation_day + isolation_period_days;
					do return_home;
				}
			}
		}

		if enable_vaccination {
			// Fractional doses are accumulated so 0.05%/day is exact over time.
			vaccination_credit <- vaccination_credit + population_size * daily_vaccination_rate;
			int doses_today <- int(vaccination_credit);
			vaccination_credit <- vaccination_credit - doses_today;
			list<person> eligible <- person where !each.vaccinated;
			doses_today <- min([doses_today, length(eligible)]);
			ask doses_today among eligible {
				vaccinated <- true;
				vaccine_variant <- vaccine_target_variant;
			}
		}

		if enable_variants {
			ask person where (each.health_state = 1) {
				if flip(daily_mutation_probability) {
					next_variant_id <- next_variant_id + 1;
					active_variant <- next_variant_id;
					infection_probability <- min([0.95, max([0.01, infection_probability * rnd(0.5, 1.5)])]);
				}
			}

			// Authorities retarget vaccination to the most prevalent active variant.
			list<int> circulating_variants <- remove_duplicates((person where (each.health_state = 1)) collect each.active_variant);
			int largest_variant_count <- -1;
			loop variant_id over: circulating_variants {
				int variant_count <- length(person where ((each.health_state = 1) and (each.active_variant = variant_id)));
				if variant_count > largest_variant_count {
					largest_variant_count <- variant_count;
					vaccine_target_variant <- variant_id;
				}
			}
		}
	}

	reflex update_peak {
		peak_infected <- max([peak_infected, infected_count]);
	}

	// GUI experiments are stopped reliably by a reflex that pauses the
	// simulation, following GAMA's official Incremental Model tutorial.
	reflex stop_simulation when: stop_when_epidemic_finished and epidemic_finished {
		do pause;
	}

}

// Lightweight procedural clouds keep the project self-contained. Their draw
// position drifts gently across the map while their agent location stays fixed.
species sky_cloud {
	float cloud_altitude <- 60#m;
	float cloud_size <- 14#m;
	float drift_speed <- 0.04#m;

	aspect three_dimensional {
		float cloud_x <- (location.x + cycle * drift_speed) mod world.shape.width;
		point cloud_center <- {cloud_x, location.y, cloud_altitude};
		rgb cloud_white <- rgb(255, 255, 255, 218);
		draw sphere(cloud_size * 0.52) at: cloud_center
			color: cloud_white lighted: true;
		draw sphere(cloud_size * 0.42)
			at: {cloud_x - cloud_size * 0.38, location.y, cloud_altitude - cloud_size * 0.08}
			color: cloud_white lighted: true;
		draw sphere(cloud_size * 0.46)
			at: {cloud_x + cloud_size * 0.38, location.y, cloud_altitude - cloud_size * 0.06}
			color: cloud_white lighted: true;
		draw sphere(cloud_size * 0.34)
			at: {cloud_x, location.y + cloud_size * 0.28, cloud_altitude - cloud_size * 0.12}
			color: rgb(235, 242, 247, 205) lighted: true;
	}
}

species road {
	// True for the single segment bridging the two city blocks across the
	// central gap. A closed connector is drawn red as a visual reminder.
	bool is_connector <- false;
	rgb road_color <- (is_connector and !enable_block_connector)
		? rgb(214, 69, 65) : rgb(108, 117, 125);

	aspect default {
		draw shape color: road_color width: 2.2#m;
	}

	// GAMA's bundled Luneray Flu and 3D GIS examples use a widened line
	// to keep the road network readable below extruded buildings.
	aspect three_dimensional {
		draw line(shape.points, 2.8#m)
			color: (is_connector and !enable_block_connector) ? rgb(214, 69, 65) : rgb(65, 72, 84)
			depth: 0.25#m;
	}
}

species building {
	bool is_school <- false;
	bool is_home <- false;
	bool is_workplace <- false;
	int community <- 0;
	float display_height <- 8#m;

	aspect default {
		rgb fill_color <- is_school ? rgb(255, 193, 7)
			: (is_workplace ? rgb(239, 71, 111) : rgb(72, 202, 228));
		draw shape color: fill_color border: rgb(82, 92, 102);
	}

	// Extruding a GIS polygon with the depth facet is the native GAMA
	// pattern used by the official Building Elevation example.
	aspect three_dimensional {
		rgb facade_color <- is_school ? rgb(255, 193, 7)
			: (is_workplace ? rgb(239, 71, 111) : rgb(72, 202, 228));
		rgb border_color <- is_school ? rgb(176, 123, 0)
			: (is_workplace ? rgb(145, 37, 77) : rgb(16, 112, 138));
		float detail_width <- min([10#m, max([4#m, shape.width * 0.42])]);
		draw shape color: facade_color border: border_color depth: display_height lighted: true;

		// Small procedural landmarks make each role readable as a cartoon
		// building while the GIS footprint remains the real building base.
		if show_cartoon_details {
			if is_home {
				// A warm pitched roof turns the cyan residential extrusion into a house.
				draw pyramid(9#m) scaled_by {1.35, 1.0, 0.42}
					at: {location.x, location.y, display_height}
					color: rgb(244, 124, 85) border: rgb(159, 68, 47) lighted: true;
				// Dark timber door, warm windows, and chimney add a lived-in scale.
				draw cube(2.2#m) scaled_by {0.75, 0.18, 1.45}
					at: {location.x, location.y - detail_width * 0.48, 1.6#m}
					color: rgb(105, 67, 48) border: rgb(69, 43, 31) lighted: true;
				draw cube(1.7#m) scaled_by {1.0, 0.14, 0.72}
					at: {location.x - 2.3#m, location.y - detail_width * 0.49, 3.0#m}
					color: rgb(255, 229, 153) border: #white lighted: false;
				draw cube(1.7#m) scaled_by {1.0, 0.14, 0.72}
					at: {location.x + 2.3#m, location.y - detail_width * 0.49, 3.0#m}
					color: rgb(255, 229, 153) border: #white lighted: false;
				draw cylinder(0.55#m, 3.4#m)
					at: {location.x + 2.5#m, location.y, display_height + 2.0#m}
					color: rgb(125, 77, 60) lighted: true;
			}
			if is_workplace {
				// A bright rooftop service block and antenna identify office towers.
				draw cube(6#m) scaled_by {1.25, 1.0, 0.45}
					at: {location.x, location.y, display_height + 1.2#m}
					color: rgb(255, 177, 199) border: rgb(145, 37, 77) lighted: true;
				draw cylinder(0.45#m, 5#m)
					at: {location.x, location.y, display_height + 4.0#m}
					color: rgb(78, 57, 86) lighted: true;
				// Reflective glass bands and a sheltered lobby modernize the offices.
				loop floor over: [0.28, 0.50, 0.72] {
					draw cube(detail_width) scaled_by {1.0, 0.08, 0.12}
						at: {location.x, location.y - detail_width * 0.52, display_height * floor}
						color: rgb(117, 205, 230) border: rgb(226, 247, 255) lighted: false;
				}
				draw cube(4.0#m) scaled_by {1.2, 0.75, 0.55}
					at: {location.x, location.y - detail_width * 0.55, 1.2#m}
					color: rgb(50, 77, 103) border: rgb(194, 229, 239) lighted: true;
				draw cube(5.5#m) scaled_by {1.25, 0.9, 0.10}
					at: {location.x, location.y - detail_width * 0.62, 3.4#m}
					color: rgb(238, 245, 249) border: rgb(145, 37, 77) lighted: true;
			}
			if is_school {
				// The school gets a central clock-tower silhouette and blue roof.
				draw cube(8#m) scaled_by {1.0, 1.0, 0.75}
					at: {location.x, location.y, display_height + 2.0#m}
					color: rgb(255, 225, 126) border: rgb(176, 123, 0) lighted: true;
				draw pyramid(9#m) scaled_by {1.0, 1.0, 0.45}
					at: {location.x, location.y, display_height + 5.0#m}
					color: rgb(56, 132, 196) border: rgb(28, 79, 121) lighted: true;
				draw sphere(1.25#m)
					at: {location.x, location.y, display_height + 4.0#m}
					color: #white border: rgb(56, 78, 98) lighted: true;
				// A blue entrance canopy, columns, and flag make the civic role clear.
				draw cube(5.0#m) scaled_by {1.5, 0.65, 0.12}
					at: {location.x, location.y - detail_width * 0.58, 4.4#m}
					color: rgb(56, 132, 196) border: rgb(28, 79, 121) lighted: true;
				draw cylinder(0.38#m, 4.2#m)
					at: {location.x - 2.4#m, location.y - detail_width * 0.52, 2.1#m}
					color: rgb(246, 246, 238) lighted: true;
				draw cylinder(0.38#m, 4.2#m)
					at: {location.x + 2.4#m, location.y - detail_width * 0.52, 2.1#m}
					color: rgb(246, 246, 238) lighted: true;
				draw cylinder(0.20#m, 11#m)
					at: {location.x + 6.0#m, location.y, 5.5#m}
					color: rgb(90, 96, 104) lighted: true;
				draw triangle(2.8#m)
					at: {location.x + 7.2#m, location.y, 10.0#m}
					color: rgb(225, 45, 55) border: rgb(142, 25, 31) lighted: false;
			}
		}

		if is_school {
			draw "SCHOOL" at: {location.x, location.y, display_height + 10.0#m}
				color: rgb(91, 61, 0) font: font("Helvetica", 13, #bold)
				anchor: #center perspective: false;
		}
		if is_workplace {
			draw "WORK" at: {location.x, location.y, display_height + 7.0#m}
				color: #white font: font("Helvetica", 8, #bold)
				anchor: #center perspective: false;
		}
	}
}

species person skills: [moving] {
	// 0 = Susceptible, 1 = Infected, 2 = Recovered.
	int health_state <- 0;
	bool is_child <- false;
	int community <- 0;
	building home;
	building day_place;
	building current_building;
	building target_building;
	point travel_target <- nil;
	bool travelling <- false;
	rgb clothing_color <- rnd_color(220);
	// Scaled walking speed for the compact 500m GIS map. With the 30-minute
	// step a commute usually completes in one or two cycles.

	int infection_start_day <- -1;
	int active_variant <- -1;
	float infection_probability <- base_infection_probability;
	list<int> immune_variants <- [];

	bool isolated <- false;
	int isolation_end_day <- -1;
	bool vaccinated <- false;
	int vaccine_variant <- 0;

	init {
		speed <- is_child ? 0.8#km/#h : rnd(0.9, 1.3)#km/#h;
	}

	action become_infected (int new_variant, float new_probability) {
		health_state <- 1;
		active_variant <- new_variant;
		infection_probability <- new_probability;
		infection_start_day <- simulation_day;
		cumulative_infections <- cumulative_infections + 1;
	}

	action start_trip (building destination_building) {
		target_building <- destination_building;
		travel_target <- any_location_in(destination_building.shape);
		current_path <- nil;
		current_building <- nil;
		travelling <- true;
	}

	action return_home {
		if current_building != home {
			do start_trip(home);
		} else {
			travelling <- false;
			target_building <- nil;
			travel_target <- nil;
		}
	}

	// With the block connector closed, a person whose daily destination lies in
	// the other community cannot cross and stays home instead of being stranded.
	reflex leave_for_day when: ((cycle mod cycles_per_day) = (departure_hour * cycles_per_hour))
		and !isolated and !travelling and (day_place != nil)
		and (enable_block_connector or (day_place.community = community)) {
		do start_trip(day_place);
	}

	reflex leave_for_home when: ((cycle mod cycles_per_day) = (return_hour * cycles_per_hour))
		and !travelling {
		do return_home;
	}

	reflex enforce_isolation when: isolated and (current_building != home) and !travelling {
		do return_home;
	}

	reflex travel when: travelling {
		do goto target: travel_target on: road_network recompute_path: false;
		if location = travel_target {
			current_building <- target_building;
			travelling <- false;
			target_building <- nil;
			travel_target <- nil;
			current_path <- nil;
		}
	}

	// Household, workplace and school contacts are evaluated twice a day.
	reflex transmit when: (health_state = 1) and !travelling and
		(((cycle mod cycles_per_day) = (work_contact_hour * cycles_per_hour)) or
		 ((cycle mod cycles_per_day) = (home_contact_hour * cycles_per_hour))) {
		list<person> contacts <- person where (
			(each != self) and !each.travelling and
			(each.current_building = current_building) and
			(each.health_state != 1) and
			((each.health_state = 0) or !(each.immune_variants contains active_variant))
		);

		// A contact cap avoids making every child meet the whole school each day.
		int contact_count <- min([max_contacts_per_event, length(contacts)]);
		list<person> met_contacts <- contact_count among contacts;
		ask met_contacts {
			float effective_probability <- myself.infection_probability;
			if vaccinated {
				// Matched vaccine: /3. A new unmatched variant is 2x less protected: /1.5.
				effective_probability <- (vaccine_variant = myself.active_variant)
					? effective_probability / 3.0
					: effective_probability / 1.5;
			}
			if flip(effective_probability) {
				do become_infected(myself.active_variant, myself.infection_probability);
			}
		}
	}

	reflex recover when: (health_state = 1) and
		(simulation_day - infection_start_day >= infectious_period_days) {
		health_state <- 2;
		if !(immune_variants contains active_variant) {
			immune_variants <- immune_variants + active_variant;
		}
		active_variant <- -1;
	}

	reflex end_isolation when: isolated and (simulation_day >= isolation_end_day) {
		isolated <- false;
	}

	aspect default {
		rgb state_color <- health_state = 0 ? rgb(46, 160, 67)
			: (health_state = 1 ? rgb(218, 54, 51) : rgb(9, 105, 218));
		float marker_size <- is_child ? 2.0#m : 2.6#m;
		draw circle(marker_size) color: state_color
			border: (isolated ? rgb(137, 87, 229) : (vaccinated ? #white : rgb(36, 41, 47)));
		if isolated {
			draw circle(marker_size + 1.1#m) color: #transparent border: rgb(137, 87, 229) width: 0.7#m;
		}
	}

	aspect three_dimensional {
		rgb state_color <- health_state = 0 ? rgb(38, 194, 129)
			: (health_state = 1 ? rgb(255, 69, 88) : rgb(53, 132, 255));
		float marker_size <- is_child ? 1.7#m : 2.2#m;
		float marker_z <- travelling ? 0.4#m
			: (current_building = nil ? 0.4#m : current_building.display_height + 0.4#m);
		point body_location <- {location.x, location.y, marker_z};
		point head_location <- {location.x, location.y, marker_z + marker_size * 0.92};
		float human_size <- is_child ? 3.0#m : 4.2#m;
		point human_location <- {location.x, location.y, marker_z + human_size * 1.35};
		point halo_location <- {location.x, location.y, marker_z + human_size * 0.78};
		point status_location <- {location.x, location.y, marker_z + human_size * 1.55};

		// Intervention halos remain visible around both realistic and cartoon
		// people, so isolation and vaccination are readable in the 3D view.
		if isolated {
			draw sphere(human_size * 0.72) at: halo_location
				color: rgb(145, 92, 246) wireframe: true lighted: true;
		} else if vaccinated {
			draw sphere(human_size * 0.62) at: halo_location
				color: rgb(255, 255, 255) wireframe: true lighted: true;
		}

		if use_realistic_people_3d {
			// Human mesh from GAMA's open Luneray Flu tutorial library.
			draw obj_file("../includes/people.obj", 90::{-1, 0, 0})
				size: human_size at: human_location rotate: heading - 90
				color: clothing_color lighted: true;
			// A small floating beacon keeps health state readable without
			// painting the entire human model red, green, or blue.
			draw sphere(is_child ? 0.38#m : 0.48#m) at: status_location
				color: state_color border: #white lighted: true;
		} else {
			draw pyramid(marker_size * 1.35) at: body_location
				color: state_color border: state_color.darker lighted: true;
			draw sphere(marker_size * 0.43) at: head_location
				color: state_color.brighter border: state_color.darker lighted: true;
		}
	}
}

experiment flu_city type: gui {
	parameter "Initial infected people" var: initially_infected min: 1 max: 100 category: "Epidemic";
	parameter "Maximum contacts per event" var: max_contacts_per_event min: 1 max: 50 category: "Epidemic";
	parameter "Workplace share of non-school buildings" var: workplace_building_fraction min: 0.05 max: 0.50 step: 0.05 category: "City roles";
	parameter "Enable families" var: enable_families category: "City roles";
	parameter "Enable school" var: enable_school category: "City roles";
	parameter "Enable block connector road" var: enable_block_connector category: "City roles";
	parameter "Show detailed 3D buildings" var: show_cartoon_details category: "City roles";
	parameter "Use realistic 3D people" var: use_realistic_people_3d category: "City roles";
	parameter "Show sky and moving clouds" var: show_sky_and_clouds category: "City roles";
	parameter "Transmission probability per contact" var: base_infection_probability min: 0.0 max: 1.0 step: 0.01 category: "Epidemic";
	parameter "Infectious period (days)" var: infectious_period_days min: 1 max: 30 category: "Epidemic";
	parameter "Stop when epidemic ends" var: stop_when_epidemic_finished category: "Epidemic";
	parameter "Enable testing and isolation" var: enable_isolation category: "Public health";
	parameter "Daily testing rate" var: daily_testing_rate min: 0.0 max: 0.20 step: 0.005 category: "Public health";
	parameter "Isolation duration (days)" var: isolation_period_days min: 1 max: 30 category: "Public health";
	parameter "Enable vaccination" var: enable_vaccination category: "Vaccination";
	parameter "Initial vaccinated coverage" var: initial_vaccination_coverage min: 0.0 max: 0.90 step: 0.10 category: "Vaccination";
	parameter "Daily vaccination rate" var: daily_vaccination_rate min: 0.0 max: 0.02 step: 0.0005 category: "Vaccination";
	parameter "Enable variants" var: enable_variants category: "Variants";
	parameter "Daily mutation probability" var: daily_mutation_probability min: 0.0 max: 0.05 step: 0.001 category: "Variants";

	output {
		monitor "Day / time" value: string(simulation_day) + " / " + string(hour_of_day) + ":" +
			(minute_of_hour = 0 ? "00" : string(minute_of_hour)) color: rgb(88, 96, 105);
		monitor "Population" value: population_size color: rgb(120, 53, 15);
		monitor "Susceptible" value: susceptible_count color: rgb(46, 160, 67);
		monitor "Infected" value: infected_count color: rgb(218, 54, 51);
		monitor "Recovered" value: recovered_count color: rgb(9, 105, 218);
		monitor "Isolated" value: isolated_count color: rgb(137, 87, 229);
		monitor "Commuting on roads" value: commuting_count color: rgb(88, 96, 105);
		monitor "Homes / workplaces / schools" value: string(length(home_buildings)) + " / " +
			string(length(workplace_buildings)) + " / " + (enable_school ? "1" : "0") color: rgb(72, 149, 239);
		monitor "Block connector" value: enable_block_connector ? "open" : "closed"
			color: enable_block_connector ? rgb(46, 160, 67) : rgb(218, 54, 51);
		monitor "Vaccinated (%)" value: vaccinated_percent color: rgb(9, 105, 218);
		monitor "Peak infected" value: peak_infected color: rgb(218, 54, 51);
		monitor "Cumulative attack rate (%)" value: attack_rate color: rgb(88, 96, 105);
		monitor "Current vaccine target" value: "Variant " + string(vaccine_target_variant) color: rgb(191, 135, 0);

		layout #split;

		display "3D Flu City" type: 3d background: rgb(112, 184, 230) antialias: true axes: false {
			light #ambient intensity: 125;
			light #default intensity: 205 direction: {0.35, 0.45, -1.0};
			graphics "Ground" refresh: false {
				draw shape color: rgb(190, 218, 184);
				// A stylized sun provides a warm visual focal point above the city.
				if show_sky_and_clouds {
					draw sphere(8#m)
						at: {world.shape.width * 0.86, world.shape.height * 0.13, 82#m}
						color: rgb(255, 221, 87) lighted: false;
				}
			}
			species sky_cloud aspect: three_dimensional visible: show_sky_and_clouds;
			species road aspect: three_dimensional refresh: false;
			species building aspect: three_dimensional refresh: false;
			species person aspect: three_dimensional;
		}

		display "Epidemic dashboard" type: 2d refresh: every(6#cycles) background: #white {
			chart "SIR population through time" type: xy style: spline
				x_label: "Simulation day" x_tick_unit: 1.0 position: {0.0, 0.0} size: {1.0, 0.62} {
				data "Susceptible" value: {simulation_time_days, susceptible_count}
					accumulate_values: true color: rgb(46, 160, 67);
				data "Infected" value: {simulation_time_days, infected_count}
					accumulate_values: true color: rgb(218, 54, 51);
				data "Recovered" value: {simulation_time_days, recovered_count}
					accumulate_values: true color: rgb(9, 105, 218);
			}
			chart "Public-health response" type: xy style: spline
				x_label: "Simulation day" x_tick_unit: 1.0 position: {0.0, 0.64} size: {1.0, 0.36} {
				data "Isolated" value: {simulation_time_days, isolated_count}
					accumulate_values: true color: rgb(137, 87, 229);
				data "Vaccinated" value: {simulation_time_days, vaccinated_count}
					accumulate_values: true color: rgb(191, 135, 0);
			}
		}
	}
}

// Extension 4: compare 10%, 50% and 90% initial vaccination coverage.
// Run this experiment in GAMA's batch view. The exhaustive method reports
// final attack rate for each coverage, repeated with controlled random seeds.
experiment vaccination_coverage_analysis type: batch repeat: 10 keep_seed: true
	until: epidemic_finished or (simulation_day >= 120) {
	parameter "Initial vaccination coverage" var: initial_vaccination_coverage among: [0.10, 0.50, 0.90];
	parameter "Testing and isolation enabled" var: enable_isolation <- true;
	parameter "Vaccination enabled" var: enable_vaccination <- true;
	parameter "Variants enabled" var: enable_variants <- true;
}

// Fast non-GUI check used to validate GIS loading and the first animated commute.
experiment automated_smoke_test type: batch repeat: 1
	until: cycle >= (11 * cycles_per_day) {
	parameter "Testing and isolation enabled" var: enable_isolation <- true;
	parameter "Vaccination enabled" var: enable_vaccination <- true;
	parameter "Variants enabled" var: enable_variants <- true;
}
