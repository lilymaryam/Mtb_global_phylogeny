import argparse

def probability(value):
    fval = float(value)
    if not 0 <= fval <= 1:
        raise argparse.ArgumentTypeError(f"{value} is not in range [0, 1]")
    return fval

def parse_arguments():
    parser = argparse.ArgumentParser(description="Subsample specific countries from the dataset.")
    parser.add_argument("--sample_country_tsv", type=str, required=True, help="Path to input tsv containing sample and 3 letter country code.")
    parser.add_argument("--output_file", type=str, required=True, help="Path to output file to be used for pruning.")
    parser.add_argument("--countries", type=str, required=True, help="Comma-separated country codes (e.g., GBR,CHN,IND)")
    parser.add_argument("--downsample", type=probability, required=True, help="Percentage of samples to remove.")
    return parser.parse_args()

def read_file(file_path):
    with open(file_path) as f:
        header = f.readline()
        countries = {}
        for line in f:
            fields = line.strip().split("\t")
            sample = fields[0]
            if len(fields) == 2:
                country_code = fields[1]
                #print(f"Sample: {sample}, Country code: {country_code}")
                if len(country_code) != 3:
                    assert False, f"Country code {country_code} is not 3 letters long."
            if len(fields) > 2:
                assert False, f"Expected 2 fields in the tsv file, but got {len(fields)} fields."
            if len(fields) == 1:
                country_code = "NA"
            if country_code not in countries:
                countries[country_code] = []
            countries[country_code].append(sample)
    return countries

def downsample_samples(countries_dict, countries, percentage, output_file):
    prune_samples = {}
    for country_code, samples in countries_dict.items():
        #print(f"Processing country code: {country_code}, Number of samples: {len(samples)}")
        if country_code in countries:
            num_samples = int(len(samples) * percentage)
            #print(f"Number of samples to prune for country", samples)
            prune_samples[country_code] = samples[:num_samples]
    
    with open(output_file, "w") as f:
        for country_code, samples in prune_samples.items():
            print(f"Country: {country_code}, Number of samples to prune: {len(samples)}")
            for sample in samples:
                f.write(sample + "\t" + country_code + "\n")


def main():
    args = parse_arguments()
    countries_dict = read_file(args.sample_country_tsv)
    country_list = args.countries.split(',')
    downsample_samples(countries_dict, country_list, args.downsample, args.output_file)

if __name__ == "__main__":
    main()


    
    
    
    