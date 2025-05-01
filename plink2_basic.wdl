version 1.0

import "https://raw.githubusercontent.com/shengqh/warp/develop/tasks/vumc_biostatistics/GcpUtils.wdl" as http_GcpUtils

workflow VUMCPlink2 {
  input {
    Array[File] source_pgen
    Array[File] source_pvar
    Array[File] source_psam
    Array[String] chromosomes

    String target_prefix
    String plink_option
    String target_gcp_folder
    String project_id

    String? parameter_file1_arg
    File? parameter_file1

    String? parameter_file2_arg
    File? parameter_file2

    String? parameter_file3_arg
    File? parameter_file3

    String docker = "shengqh/plink_1.9_2.0:20250304"
    Int? memory_size = 10
  }

  scatter (idx in range(length(chromosomes))) {
    String chromosome = chromosomes[idx]
    File my_source_pgen = source_pgen[idx]
    File my_source_pvar = source_pvar[idx]
    File my_source_psam = source_psam[idx]

    String my_target_prefix = target_prefix + "_" + chromosome

    call Plink2 {
      input:
        source_pgen = my_source_pgen,
        source_pvar = my_source_pvar,
        source_psam = my_source_psam,
        plink_option = plink_option,
        out_string = my_target_prefix,
        parameter_file1_arg = parameter_file1_arg,
        parameter_file1 = parameter_file1,
        parameter_file2_arg = parameter_file2_arg,
        parameter_file2 = parameter_file2,
        parameter_file3_arg = parameter_file3_arg,
        parameter_file3 = parameter_file3,
        docker = docker,
        memory_size = memory_size
    }
  }

  call MergePgenFiles {
    input:
      pgen_files = Plink2.output_pgen,
      pvar_files = Plink2.output_pvar,
      psam_files = Plink2.output_psam,
      output_prefix = target_prefix
  }

  call http_GcpUtils.MoveOrCopyThreeFiles as CopyFiles_two {
    input:
      source_file1 = MergePgenFiles.output_pgen,
      source_file2 = MergePgenFiles.output_pvar,
      source_file3 = MergePgenFiles.output_psam,
      is_move_file = false,
      project_id = project_id,
      target_gcp_folder = target_gcp_folder
  }

  output {
    Array[File] ancestry_outputs = [
      CopyFiles_two.output_file1,
      CopyFiles_two.output_file2,
      CopyFiles_two.output_file3
    ]
  }
}

task Plink2 {
  input {
    File source_pgen
    File source_pvar
    File source_psam

    String plink_option
    String out_string

    String? parameter_file1_arg
    File? parameter_file1

    String? parameter_file2_arg
    File? parameter_file2

    String? parameter_file3_arg
    File? parameter_file3

    String docker = "shengqh/plink_1.9_2.0:20250304"
    Int? memory_size = 10
  }

  Int disk_size = ceil(size(source_pgen, "GB") * 2) + 2

  command <<<

    plink2 \
      --pgen ~{source_pgen} \
      --pvar ~{source_pvar} \
      --psam ~{source_psam} \
      ~{parameter_file1_arg + " " + parameter_file1} \
      ~{parameter_file2_arg + " " + parameter_file2} \
      ~{parameter_file3_arg + " " + parameter_file3} \
      ~{plink_option} \
      --out ~{out_string}

  >>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_size + " GiB"
  }

  output {
    File output_pgen = out_string + ".pgen"
    File output_pvar = out_string + ".pvar"
    File output_psam = out_string + ".psam"
  }
}

task MergePgenFiles {
  input {
    Array[File] pgen_files
    Array[File] pvar_files
    Array[File] psam_files

    String output_prefix

    Int memory_gb = 20
    Int cpu = 8

    String docker = "shengqh/plink_1.9_2.0:20250304"
  }

  Int disk_size = ceil((size(pgen_files, "GB") + size(pvar_files, "GB") + size(psam_files, "GB")) * 3) + 20

  String target_pgen = output_prefix + ".pgen"
  String target_pvar = output_prefix + ".pvar"
  String target_psam = output_prefix + ".psam"

  String merged_pgen = output_prefix + "-merge.pgen"
  String merged_pvar = output_prefix + "-merge.pvar"
  String merged_psam = output_prefix + "-merge.psam"

  command <<<

    cat ~{write_lines(pgen_files)} > pgen.list
    cat ~{write_lines(pvar_files)} > pvar.list
    cat ~{write_lines(psam_files)} > psam.list

    paste pgen.list pvar.list psam.list > merge.list

    plink2 --pmerge-list merge.list --make-pgen --out ~{output_prefix} --threads ~{cpu}

    rm -f ~{target_pgen} ~{target_pvar} ~{target_psam}

    mv ~{merged_pgen} ~{target_pgen}
    mv ~{merged_pvar} ~{target_pvar}
    mv ~{merged_psam} ~{target_psam}

    grep -v "^#" ~{target_psam} | wc -l | cut -d ' ' -f 1 > num_samples.txt
    grep -v "^#" ~{target_pvar} | wc -l | cut -d ' ' -f 1 > num_variants.txt

  >>>

  runtime {
    cpu: cpu
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }

  output {
    File output_pgen = target_pgen
    File output_pvar = target_pvar
    File output_psam = target_psam

    Int num_samples = read_int("num_samples.txt")
    Int num_variants = read_int("num_variants.txt")
  }
}