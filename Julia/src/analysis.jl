function analyze(Results::DataFrame)
    mean_and_std = 
    combine(groupby(Results, :temperature), 
        :average_energy => mean, :average_energy => std,
        :magnetisation  => mean, :magnetisation  => std,
        :molecules      => mean, :molecules      => std
        );
    
    return mean_and_std;
end

