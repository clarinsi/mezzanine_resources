from pathlib import Path
dirs_to_process = [
    "iriss",
    "SST",
    "SPOG"
]

files_to_process = []
for i in dirs_to_process:
    files_to_process += list(Path(i).glob("**/*.xml"))
wd = "snakedir"
rule PREP_WORKING_DIR:
    output: [directory(f"{wd}/{i}") for i in dirs_to_process]
    shell:
        f"mkdir -p {wd} " + " ".join([f"{wd}/{i}" for i in dirs_to_process])

rule DO_ONE:
    input: "{dir}/{file}.xml" , [f"{wd}/{i}" for i in dirs_to_process]
    output: f"{wd}/{{dir}}/{{file}}.csv"
    conda: "mezzanine.yml"
    script:
        "doone.py"

rule GATHER_RESULTS:
    input: [f"{wd}/{i.with_suffix('.csv')}" for i in files_to_process]
    output: "calculation.done"
    shell:
        "touch calculation.done"

rule AGGREGATE_SUBCORPUS:
    input: f"{wd}/{{dir}}/", "calculation.done"
    output:
        csv=f"{{wd}}/{{dir}}.csv",
        html=f"{{wd}}/{{dir}}.html",
    params:
        template="template.html",
    conda: "mezzanine.yml"
    script:
        "gather.py"

rule GATHER_AGGREGATES:
    default_target: True
    input: [f"{wd}/{i}.csv" for i in dirs_to_process]


rule clean:
    shell:
        f"rm -rf {wd} calculation.done"