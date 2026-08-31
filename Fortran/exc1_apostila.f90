PROGRAM exc1_apostila  
! Purpose:  Thermodynamics exercise 1 from the apostila! A program about bissection method!
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
    REAL(dp) :: a, b , c
    REAL(dp) :: media, erro
    INTEGER :: i

! 1. Prompt the user for the input
    


! 2. Calculations

    ! 2.2 Initialize variables

    i = 0 
    a = 300.0d0
    b = 500.0d0   
    erro = ABS(b - a)

    ! Evaluate the function at the endpoints to check for a sign change
    IF (f(a) * f(b) > 0.0d0) THEN
        WRITE(*,*) "Erro: O intervalo [a, b] nao contem troca de sinal em f(T)."
        STOP
    END IF

    DO WHILE (erro >= tol)
        
        c = 0.5d0 * (a + b)
        
        IF (f(a) * f(c) >= 0.0d0) THEN
            a = c
        ELSE
            b = c
        END IF
        
        erro = ABS(b - a)
        i = i + 1
    END DO

    ! Convergency Analysis 
    IF (erro < tol) THEN
        WRITE(*,*) "The method converged successfully."
    ELSE
        WRITE(*,*) "The method did not converge within the specified tolerance."
    END IF

    ! Calculate the midpoint of the final interval
    media = 0.5d0 * (a + b)    

! 3. Write out the result

    WRITE(*,*) "The root is: ", media
    WRITE(*,*) "Number of iterations: ", i
    WRITE(*,*) "Final interval: [", a, ", ", b, "]"
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

END PROGRAM exc1_apostila