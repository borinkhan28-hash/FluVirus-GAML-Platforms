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
	// Ten-minute cycles keep road travel visible instead of teleporting agents.
	float step <- 10#mn;
	int cycles_per_hour <- 6;
	int cycles_per_day <- 144;

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
	float initial_vaccination_coverage <- 0.10;
	float daily_vaccination_rate <- 0.0005;
	float daily_mutation_probability <- 0.001;

	// Runtime state and indicators.
	int simulation_day -> int(cycle / cycles_per_day);
	int hour_of_day -> int((cycle mod cycles_per_day) / cycles_per_hour);
	int minute_of_hour -> int((cycle mod cycles_per_hour) * 10);
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
		road_network <- as_edge_graph(road);

		// Building-role logic follows the separated-role reference model:
		// one school, a configurable workplace group, and family homes.
		school <- building with_max_of (each.shape.area);
		ask school {
			is_school <- true;
			display_height <- 20#m;
		}
		list<building> available_buildings <- building where (each != school);
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

		// Every residential building receives one family of 3-6 people
		// with 0-2 children. Adults work; children attend the school.
		ask home_buildings {
			int family_size <- rnd(3, 6);
			int child_count <- rnd(0, min([2, family_size - 1]));

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

}

species road {
	aspect default {
		draw shape color: rgb(108, 117, 125) width: 2.2#m;
	}

	// GAMA's bundled Luneray Flu and 3D GIS examples use a widened line
	// to keep the road network readable below extruded buildings.
	aspect three_dimensional {
		draw line(shape.points, 2.8#m) color: rgb(65, 72, 84) depth: 0.25#m;
	}
}

species building {
	bool is_school <- false;
	bool is_home <- false;
	bool is_workplace <- false;
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
		draw shape color: facade_color border: border_color depth: display_height lighted: true;

		// Small procedural landmarks make each role readable as a cartoon
		// building while the GIS footprint remains the real building base.
		if show_cartoon_details {
			if is_home {
				// A warm pitched roof turns the cyan residential extrusion into a house.
				draw pyramid(9#m) scaled_by {1.35, 1.0, 0.42}
					at: {location.x, location.y, display_height}
					color: rgb(244, 124, 85) border: rgb(159, 68, 47) lighted: true;
			}
			if is_workplace {
				// A bright rooftop service block and antenna identify office towers.
				draw cube(6#m) scaled_by {1.25, 1.0, 0.45}
					at: {location.x, location.y, display_height + 1.2#m}
					color: rgb(255, 177, 199) border: rgb(145, 37, 77) lighted: true;
				draw cylinder(0.45#m, 5#m)
					at: {location.x, location.y, display_height + 4.0#m}
					color: rgb(78, 57, 86) lighted: true;
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
	building home;
	building day_place;
	building current_building;
	building target_building;
	point travel_target <- nil;
	bool travelling <- false;
	rgb clothing_color <- rnd_color(220);
	// Scaled walking speed gives several visible frames on the compact 500m GIS map.

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

	reflex leave_for_day when: ((cycle mod cycles_per_day) = (departure_hour * cycles_per_hour))
		and !isolated and !travelling {
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

		// A pyramid torso plus spherical head gives agents a readable cartoon
		// person silhouette. The halo preserves intervention information.
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
	parameter "Show cartoon building details" var: show_cartoon_details category: "City roles";
	parameter "Use realistic 3D people" var: use_realistic_people_3d category: "City roles";
	parameter "Transmission probability per contact" var: base_infection_probability min: 0.0 max: 1.0 step: 0.01 category: "Epidemic";
	parameter "Infectious period (days)" var: infectious_period_days min: 1 max: 30 category: "Epidemic";
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
		monitor "Susceptible" value: susceptible_count color: rgb(46, 160, 67);
		monitor "Infected" value: infected_count color: rgb(218, 54, 51);
		monitor "Recovered" value: recovered_count color: rgb(9, 105, 218);
		monitor "Isolated" value: isolated_count color: rgb(137, 87, 229);
		monitor "Commuting on roads" value: commuting_count color: rgb(88, 96, 105);
		monitor "Homes / workplaces / schools" value: string(length(home_buildings)) + " / " +
			string(length(workplace_buildings)) + " / 1" color: rgb(72, 149, 239);
		monitor "Vaccinated (%)" value: vaccinated_percent color: rgb(9, 105, 218);
		monitor "Peak infected" value: peak_infected color: rgb(218, 54, 51);
		monitor "Cumulative attack rate (%)" value: attack_rate color: rgb(88, 96, 105);
		monitor "Current vaccine target" value: "Variant " + string(vaccine_target_variant) color: rgb(191, 135, 0);

		layout #split;

		display "3D Flu City" type: 3d background: rgb(232, 244, 248) antialias: true {
			light #ambient intensity: 110;
			light #default intensity: 190 direction: {0.5, 0.5, -1.0};
			graphics "Ground" refresh: false {
				draw shape color: rgb(221, 235, 226);
			}
			species road aspect: three_dimensional refresh: false;
			species building aspect: three_dimensional refresh: false;
			species person aspect: three_dimensional;
		}

		display "Epidemic dashboard" type: 2d refresh: every(6#cycles) background: #white {
			chart "SIR population through time" type: series style: spline position: {0.0, 0.0} size: {1.0, 0.62} {
				data "Susceptible" value: susceptible_count color: rgb(46, 160, 67);
				data "Infected" value: infected_count color: rgb(218, 54, 51);
				data "Recovered" value: recovered_count color: rgb(9, 105, 218);
			}
			chart "Public-health response" type: series style: spline position: {0.0, 0.64} size: {1.0, 0.36} {
				data "Isolated" value: isolated_count color: rgb(137, 87, 229);
				data "Vaccinated" value: vaccinated_count color: rgb(191, 135, 0);
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
