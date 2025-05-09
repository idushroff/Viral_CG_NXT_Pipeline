#!/usr/bin/env nextflow

/*testing - use echo to print 'Hello World!' to a file*/ 

process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output: 
        path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
    
}

//Pipeline parameters

params.greeting = "Hola"

//Pipeline workflow

workflow {

    greeting_ch = Channel.of('Hello Channels', 'Hello', 'Hola')
    //emit a greeting
    sayHello(greeting_ch)
}