capture program drop tt8
program tt8,rclass

    version 13

    syntax , Matt(string)
    tempname A B

    mat `A' = (1,2,3,4)

    mat `B' = `A'+`matt'

    //mat B = `B'
	mat list `B'
 end
