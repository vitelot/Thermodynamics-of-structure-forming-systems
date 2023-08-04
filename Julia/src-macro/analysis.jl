function analyze(Results::DataFrame)
    mean_and_std = 
    combine(groupby(Results, :temperature), 
        :average_energy => mean, :average_energy => std,
        :magnetisation  => mean, :magnetisation  => std,
        :molecules      => mean, :molecules      => std
        );
    
    return sort(mean_and_std, :temperature);
end

function plotall(data::DataFrame)

    p1 = plot(data.temperature, data.average_energy_mean, yerr=data.average_energy_std,
              ylabel="energy");
    p2 = plot(data.temperature, abs.(data.magnetisation_mean), yerr=data.magnetisation_std,
              ylabel = "magnetisation");
    p3 = plot(data.temperature, data.molecules_mean, yerr=data.molecules_std,
              ylabel = "molecules");

    plot(p1,p2,p3, label="", xguide="temperature")
end
