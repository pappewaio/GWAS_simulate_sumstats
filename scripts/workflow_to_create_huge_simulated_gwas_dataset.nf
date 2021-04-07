nextflow.enable.dsl=2

process foo {
    output:
      path 'foo.txt'
    script:
      """
      echo "hej" > foo.txt
      """
}

 process bar {
    input:
      path x
    output:
      path 'bar.txt'
    script:
      """
      cat $x > bar.txt
      """
}

workflow {
    foo()
    bar(foo.out)
}


