MODULE modul
  IMPLICIT NONE
  INTEGER :: n_p                 ! Nombre de particules
  REAL    :: m_part , m_brown     ! Masses relatives
  REAL    :: l_dom , r_brown       ! Tailles relatives domaine/macroparticule
  REAL    :: v_c                        ! Vitesse caracteristique des particules
  REAL    :: t_fin                     ! Temps final
  INTEGER, PARAMETER :: n_aff = 40         ! Taille de l'affichage
  LOGICAL :: affichage            ! affichage on/off
TYPE particule
    REAL :: x, y    ! position
    REAL :: vx, vy  ! vitesse
    REAL :: m       ! masse
  END TYPE particule
 
  TYPE macroparticule
    REAL :: x, y    ! position
    REAL :: vx, vy  ! vitesse
    REAL :: m, r    ! masse et taille (rayon)
  END TYPE macroparticule
 
  CONTAINS
  LOGICAL FUNCTION isin(p, b)
  ! Cette fonction retourne VRAI si la particule p est situee
  ! a l'interieur de la macroparticule b, FAUX sinon
    TYPE(particule), INTENT(IN)      :: p
    TYPE(macroparticule), INTENT(IN) :: b

    IF ((p%x-b%x)**2 + (p%y-b%y)**2 < b%r**2) THEN
      isin = .TRUE.
    ELSE
      isin = .FALSE.
    END IF
  END FUNCTION isin
   
  SUBROUTINE affiche(part, b)
  ! Ce sous-programme affiche le domaine carre, les particules
  ! du tableau part et la macroparticule b
    TYPE(particule), DIMENSION(:), INTENT(INOUT) :: part
    TYPE(macroparticule), INTENT(IN)             :: b
    CHARACTER, DIMENSION(-1:n_aff+1, -1:n_aff+1) :: tabc
    INTEGER                                      :: ip, i, j
   
    ! Effacage de l'ecran, initialisation tableau et bordures
    CALL SYSTEM('clear')
    tabc = ' '
    tabc(-1:n_aff+1:n_aff+2,:)="*"
    tabc(:,-1:n_aff+1:n_aff+2)="*"

    ! placement des particules
    DO ip = 1, n_p
      j = NINT(part(ip)%x/l_dom*n_aff)
      i = n_aff-NINT(part(ip)%y/l_dom*n_aff)
      tabc(i,j) = ACHAR(ip+64)!"o"
    END DO
    ! placement de la macroparticule
    DO i = 0, n_aff
      DO j = 0, n_aff
        IF ((l_dom - l_dom*i/n_aff - b%y)**2 + (l_dom*j/n_aff - b%x)**2 < b%r**2) tabc(i,j) = "@"
      END DO
    END DO
    ! affichage proprement dit
    DO i = -1, n_aff+1
      PRINT "(1000A2)", tabc(i, :)
    END DO
    CALL SYSTEM('sleep 0.1')
  END SUBROUTINE affiche 

SUBROUTINE readParam()

  OPEN(UNIT=343, FILE='param.dat', ACTION="READ", STATUS="OLD")
  READ(343,*) n_p, m_part, m_brown,l_dom,r_brown,v_c,t_fin,affichage  
  close(343)

END SUBROUTINE readParam


SUBROUTINE Initialisations(part, brown,n_t,dt,vmax)

    TYPE(particule), DIMENSION(:), ALLOCATABLE :: part
    TYPE(macroparticule)                       :: brown
    INTEGER                                    :: i, ok
    INTEGER                                    :: n_t,ip
    REAL                                       :: vmax, dt,u1,u2,vpx,vpy
   ! Initialisation de la macroparticule
  brown = macroparticule(0.5*l_dom, 0.5*l_dom, 0., 0., m_brown, r_brown)
    ! Initialisation des microparticules
  ALLOCATE(part(n_p), STAT=ok)
  IF (ok/=0) STOP "Echec allocation tableau particules !"
       ! positions
  CALL RANDOM_NUMBER(part%x)
  CALL RANDOM_NUMBER(part%y)
  part%x = part%x*l_dom
  part%y = part%y*l_dom
        ! vitesses

  CALL RANDOM_NUMBER(u1)
  CALL RANDOM_NUMBER(u2)
  vpx=v_c*SQRT(-2*log10(u1))*cos(2*3.14*u2)
  vpy=v_c*SQRT(-2*log10(u1))*sin(2*3.14*u2)
  Do while(SQRT(vpx**2+vpy**2)>=3*v_c  )
  CALL RANDOM_NUMBER(u1)
  CALL RANDOM_NUMBER(u2)
  vpx=v_c*SQRT(-2*log10(u1))*cos(2*3.14*u2)
  vpy=v_c*SQRT(-2*log10(u1))*sin(2*3.14*u2)
  END DO
  part%vx = vpx
  part%vy = vpy
        ! masses
  part%m = m_part
        ! rectification de la position si trop proche de brown
  DO i = 1, n_p
    DO WHILE (isin(part(i),brown))
      CALL RANDOM_NUMBER(part(i)%x)
      CALL RANDOM_NUMBER(part(i)%y)
      part(i)%x = part(i)%x*l_dom
      part(i)%y = part(i)%y*l_dom
    END DO
  END DO
    ! Choix du pas de temps
  vmax = MAXVAL(SQRT(part%vx**2+part%vy**2)) ! Vitesse max
  dt = 0.1*r_brown/vmax
  n_t = CEILING(t_fin/dt) ! Nombre de pas de temps



END SUBROUTINE Initialisations



SUBROUTINE OutputFile(trajectoireX,trajectoireY,n_t)
  INTEGER                                    :: n_t, it
  REAL, DIMENSION(:), ALLOCATABLE :: trajectoireX,trajectoireY
  OPEN(UNIT=20, FILE='output.dat', ACTION="WRITE", STATUS="replace")
  DO it = 1, n_t   
  WRITE(20,*) trajectoireX(it),trajectoireY(it)
  END DO
  close(20)
END SUBROUTINE OutputFile

SUBROUTINE  Eloignement(trajectoireX,trajectoireY,n_t,dt)
  INTEGER                                    :: n_t, it
  REAL, DIMENSION(:), ALLOCATABLE :: trajectoireX,trajectoireY
  REAL				  :: rlt,dt 

  OPEN(UNIT=20, FILE='Eloignement.dat', ACTION="WRITE", STATUS="replace")
  DO it = 1, n_t   
  rlt=rlt+ sqrt(trajectoireX(it)**2+trajectoireY(it)**2)**2
  WRITE(20,*) it*dt,rlt

  END DO
  close(20)
END SUBROUTINE Eloignement



SUBROUTINE Simulation(part, brown,n_t,dt,vmax,trajectoireX,trajectoireY)


    TYPE(particule), DIMENSION(:), ALLOCATABLE :: part
    TYPE(macroparticule)                       :: brown
    INTEGER                                    :: it, i,ok
    INTEGER                                    :: n_t,ip
    REAL                                       :: vmax, dt
    REAL :: cx, cy,Thetac, tmp
    REAL :: vpx, vpy
    REAL :: vbx, vby
    REAL :: vcpx, vcpy  !  les vitesses apres choc
    REAL :: vcbx, vcby  !  les vitesses apres choc
    REAL, DIMENSION(:), ALLOCATABLE :: trajectoireX,trajectoireY


    ALLOCATE (trajectoireX(0:n_t), stat=ok)
    IF (ok /= 0) STOP "Pblm allocation tab"

    ALLOCATE (trajectoireY(0:n_t), stat=ok)
    IF (ok /= 0) STOP "Pblm allocation tab"

  ! Dynamique moleculaire

    ! Boucle en temps
  DO it = 1, n_t
      ! Progression des particules

    part%x = part%x + part%vx*dt
    part%y = part%y + part%vy*dt
    brown%x = brown%x + brown%vx*dt
    brown%y = brown%y + brown%vy*dt 

    ! Traitement des chocs contre la paroi
   
    WHERE (part%x < 0 .OR. part%x > l_dom) part%vx = -part%vx
    WHERE (part%y < 0 .OR. part%y > l_dom) part%vy = -part%vy
    IF (brown%x - brown%r < 0 .OR. brown%x + brown%r > l_dom) brown%vx = -brown%vx
    IF (brown%y - brown%r < 0 .OR. brown%y + brown%r > l_dom) brown%vy = -brown%vy


    !Chocs microparticule/macroparticule
    DO ip = 1, n_p
    tmp= (part(ip)%vy-brown%vy)*(part(ip)%y-brown%y) +(part(ip)%vx-brown%vx)*( part(ip)%x-brown%x )
    IF ( tmp <0 .AND. isin(part(ip), brown))  THEN  
 
    cx=part(ip)%x
    cy=part(ip)%y
    Thetac  = atan2(brown%x-cx,brown%y-cy)
    vpx= cos (Thetac) *part(ip)%vx + sin(Thetac) * part(ip)%vy
    vpy= -sin (Thetac) *part(ip)%vx + cos(Thetac) * part(ip)%vy
    vbx= cos (Thetac) *brown%vx + sin(Thetac) * brown%vy
    vby= -sin (Thetac) *brown%vx + cos(Thetac) * brown%vy
    vcpx=1./(part(ip)%m + brown%m) *( (part(ip)%m - brown%m)*vpx + 2*brown%m* vbx )
    vcbx=1./(part(ip)%m + brown%m) *( 2*part(ip)%m* vpx + ( brown%m - part(ip)%m  )*vbx)
    vcpy=vpy
    vcby=vby
    part(ip)%vx= cos (Thetac) * vcpx -sin(Thetac)*vcpy  
    part(ip)%vy= sin (Thetac) * vcpx +cos(Thetac)*vcpy
    brown%vx= cos (Thetac) * vcbx -sin(Thetac)*vcby
    brown%vy= sin (Thetac) * vcbx +cos(Thetac)*vcby

    END IF

    END DO  
    trajectoireX(it)=brown%x
    trajectoireY(it) =brown%y
      ! Affichage
    IF (affichage) CALL affiche(part, brown)
  END DO


END SUBROUTINE Simulation


       
END MODULE modul







PROGRAM brownian
! Ce programme simule l'evolution temporelle d'une particule brownienne
! (brown) sous l'effet des chocs de microparticules sans taille (tableau part)
  USE modul
  IMPLICIT NONE
  TYPE(particule), DIMENSION(:), ALLOCATABLE :: part
  TYPE(macroparticule)                       :: brown
  REAL, DIMENSION(:), ALLOCATABLE            :: trajectoireX,trajectoireY,Nbiterations
  INTEGER                                    :: n_t,ip
  REAL                                       :: vmax, dt
 
  CALL readParam()
  CALL Initialisations(part, brown,n_t,dt,vmax)
  CALL Simulation(part, brown,n_t,dt,vmax,trajectoireX,trajectoireY)
  CALL OutputFile(trajectoireX,trajectoireY,n_t)

  CALL Eloignement(trajectoireX,trajectoireY,n_t,dt)

    ! Fin boucle en temps
  DEALLOCATE(part)
END PROGRAM brownian  
