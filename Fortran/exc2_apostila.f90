PROGRAM exc2_apostila  
! Purpose:  Thermodynamics exercise 2 from the apostila! A program about Newton-Raphson method!
!           

! Record of revisions:
! Date              Programmer              Description of change
! ====              ==========              =====================
! 20/08/2026 --    Daniel Campiotti            Original code
!
    USE iso_fortran_env, ONLY: dp => real64
    IMPLICIT NONE 

! Data dictionary: declare variable types, definitions, & units
    REAL(dp), PARAMETER :: tol = 1.0d-10
    REAL(dp) :: a, y, derivative
    REAL(dp) :: erro
    INTEGER :: i

! 1. Prompt the user for the input
    


! 2. Calculations

    ! 2.2 Initialize variables

    i = 0
    a = 500.0d0
    y = f(a)
    derivative = df(a)
    erro = ABS(y)

    DO WHILE (erro >= tol)

        derivative = df(a)
        a = a - (y / derivative)
        y = f(a)                    ! important to recalculate y after updating a so the error is calculated correctly
        erro = ABS(y)
        i = i + 1

    END DO
    
    ! Convergency Analysis 
    IF (erro < tol) THEN
        WRITE(*,*) "The method converged successfully."
    ELSE
        WRITE(*,*) "The method did not converge within the specified tolerance."
    END IF
 
! 3. Write out the result

    WRITE(*,*) "Number of iterations: ", i
    write(*,*) "Final value of a: ", a
    WRITE(*,*) "Final value of f(a): ", y
    write(*,*) "Final value of derivative: ", derivative
    WRITE(*,*) "Final error: ", erro


    ! Functions definitions
    CONTAINS
        FUNCTION f(T_2)
            
            IMPLICIT NONE 

            REAL(dp), INTENT(IN) :: T_2
            REAL(dp) :: f
            REAL(dp) :: R, p_1, p_2, T_1
            REAL(dp) :: a_0, a_1, a_2, a_3, a_4

            ! Contants
            R = 8.314d0     
            T_1 = 500.0d0
            p_1 = 10.0d0
            p_2 = 1.0d0

            ! Coefficients
            a_0 = 3.259d0
            a_1 = 1.356d-3
            a_2 = 1.502d-5
            a_3 = -2.374d-8
            a_4 = 1.056d-11

            f = R * ((a_0 * LOG(T_2) + a_1  * T_2 + (a_2/2.0e0) * T_2**2 + (a_3/3.0e0) * T_2**3 + (a_4 / 4.0e0) * T_2**4)&
             - (a_0 * LOG(T_1) + a_1 * T_1 + (a_2/2.0e0) * T_1**2 + (a_3/3.0e0) * T_1**3 + (a_4 / 4.0e0) * T_1**4))&
             - R * LOG(p_2 / p_1)

        END FUNCTION f
        
        FUNCTION df(T_2)
            
            IMPLICIT NONE 

            REAL(dp), INTENT(IN) :: T_2
            REAL(dp) :: df
            REAL(dp) :: R, a_0, a_1, a_2, a_3, a_4

            ! Contants
            R = 8.314d0     

            ! Coefficients
            a_0 = 3.259d0
            a_1 = 1.356d-3
            a_2 = 1.502d-5
            a_3 = -2.374d-8
            a_4 = 1.056d-11

            df = R * (a_0 / T_2 + a_1 + a_2 * T_2 + a_3 * T_2**2 + a_4 * T_2**3)

        END FUNCTION df

END PROGRAM exc2_apostila