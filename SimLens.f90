    !Simple program demonstrating how to generate a simulated lensed map
    !AL, Feb 2004; Updated Oct 2007
    program SimLensCMB
    use HealpixObj
    use HealpixVis
    use Random
    use spinalm_tools
    use IniFile
    use AMLUtils
    implicit none
    Type(HealpixInfo)  :: H
    Type(HealpixMap)   :: M, GradPhi
    Type(HealpixPower) :: P
    Type(HealpixAlm)   :: A
    integer            :: nside, lmax
    integer(I_NPIX)    :: npix
    character(LEN=1024)  :: w8name = '../Healpix_2.00/data/'
    character(LEN=1024)  :: file_stem, cls_file, out_file_root, cls_lensed_file
    character(LEN=1024) :: healpixloc
    integer, parameter :: lens_interp =1, lens_exact = 2
    integer :: lens_method = lens_interp
    integer :: mpi_division_method = division_equalrows
    integer ::  interp_method,  rand_seed
    logical :: err, want_pol
    real :: interp_factor
    integer status
    !LP
    integer npol_in
    Type(HealpixMap) :: unlens_m
    Type(HealpixAlm) :: phi_alm, a_temp
    character(LEN=1024) :: phi_file, primary_file, unlensed_cls_file, unlensed_map_file, lensed_cls_file, lensed_map_file
    logical :: unlensed_map_file_exists, lensed_map_file_exists, test_unlensed_exists
    !/LP
#ifdef MPIPIX
    integer i
    
    call mpi_init(i)
#endif
    Ini_Fail_On_Not_Found = .true.
    call Ini_Open(GetParam(1), 3,err)
    if (err) then
#ifdef MPIPIX
        call mpi_finalize(i)
#endif
        stop 'No ini'
    end if
    nside  = Ini_Read_Int('nside')
    npix = nside2npix(nside)

    lmax   = Ini_Read_Int('lmax')  
    cls_file = Ini_Read_String('cls_file')
    out_file_root = Ini_Read_String('out_file_root')

    lens_method = Ini_Read_Int('lens_method')
    want_pol = Ini_Read_Logical('want_pol')
    rand_seed = Ini_Read_Int('rand_seed')

    interp_method = Ini_read_int('interp_method')

    Ini_Fail_On_Not_Found = .false.

    !***read primary filename
    primary_file = Ini_Read_String('primary_file')
    !***read phi filename
    phi_file = Ini_Read_String('phi_file')

    w8name = Ini_Read_String('w8dir')
    interp_factor=0
    if (lens_method == lens_interp) interp_factor = Ini_Read_Real('interp_factor',3.)
#ifdef MPIPIX
    mpi_division_method = Ini_Read_Int('mpi_division_method',division_balanced);
#endif 

    call Ini_Close

    file_stem =  trim(out_file_root)//'_lmax'//trim(IntToStr(lmax))//'_nside'//trim(IntTOStr(nside))// &
    '_interp'//trim(RealToStr(interp_factor,3))//'_method'//trim(IntToStr(interp_method))//'_'

    if (want_pol) file_stem=trim(file_stem)//'pol_'
    file_stem = trim(file_stem)//trim(IntToStr(lens_method)) 
    
    call SetIdlePriority()

    if (w8name=='') then
        call get_environment_variable('HEALPIX', healpixloc, status=status)
        if (status==0) then
            w8name = trim(healpixloc)//'/data/'
        end if
    end if

    if (w8name=='') then
        write (*,*) 'Warning: using unit weights as no w8dir found'
        call HealpixInit(H,nside, lmax,.true., w8dir='', method= mpi_division_method) 
    else
        call HealpixInit(H,nside, lmax,.true., w8dir=w8name,method=mpi_division_method) 
    end if 

    !***file names to save maps/powerspectra to
    unlensed_map_file = trim(file_stem)//"_unlensed.fits"
    lensed_map_file = trim(file_stem)//"_lensed.fits"
    unlensed_cls_file = trim(file_stem)//"_unlensed_power.dat"
    lensed_cls_file = trim(file_stem)//"_lensed_power.dat"

    if (want_pol) npol_in=3

    if (H%MpiID ==0) then !if we are main thread
        !All but main thread stay in HealpixInit

        call HealpixPower_nullify(P)
        call HealpixAlm_nullify(A)
        call HealpixMap_nullify(GradPhi)
        call HealpixMap_nullify(M)
        call HealpixAlm_nullify(phi_alm)
        call HealpixAlm_nullify(a_temp)
        call HealpixMap_nullify(unlens_m)
        call HealpixPower_ReadFromTextFile(P,cls_file,lmax,pol=.true.,dolens = .true.)
        !Reads in unlensed C_l text files as produced by CAMB (or CMBFAST if you aren't doing lensing)

        !***read in primary
        write(*,*) "reading primary..."
        !call HealpixMap_Read(M, "test_primary.fits"
        !call HealpixMap2Alm(H, M, A, lmax, dopol = want_pol) 
        call HealpixAlm_Read(A, primary_file, lmax=lmax, npol_in=npol_in)
        !write(*,*) "npol: ", A%npol
        write(*,*) "primary read successful"
        call HealpixAlm2Map(H, A, M, npix)
        
        !***read in phi map (alm)
        write(*,*) "reading phi map..."
        call HealpixAlm_Read(phi_alm, phi_file, lmax=lmax)
        write(*,*) "phi read successful"
        !***convert to gradient map                                                                
        call HealpixAlm2GradientMap(H, phi_alm, GradPhi,npix,'TEB')
        
        !call HealpixAlm_Sim(A, P, rand_seed,HasPhi=.true., dopol = want_pol)
        call HealpixAlm2Power(A,P)
        call HealpixPower_Write(P,unlensed_cls_file)
        !call HealpixAlm2GradientMap(H,A, GradPhi,npix,'PHI')                                
        
        !***write unlensed map (as alm) to file -- won't need when supplying own primary
        !inquire(file=unlensed_map_file, exist=unlensed_map_file_exists)
        !if (unlensed_map_file_exists) call DeleteFile(unlensed_map_file)
        !call HealpixAlm_Write(A, unlensed_map_file)
        
        !***TEST - output unlensed as map instead of alm
        !NOTE -- can be saved, but causes lensing part to error
        !call HealpixAlm2Map(H, A, unlens_m, npix)
        !call HealpixMap_Write(unlens_m, unlensed_map_file, overwrite=.true.)

        if (lens_method == lens_exact) then
            call HealpixExactLensedMap_GradPhi(H,A,GradPhi,M) !(H,A,GradPhi,M)
        else if (lens_method == lens_interp) then
           write(*,*) "interp start"
           call HealpixInterpLensedMap_GradPhi(H,A,GradPhi, M, interp_factor, interp_method) 
           write(*,*) "interp end"
        else
            stop 'unknown lens_method'
        end if
               
        call HealpixMap2Alm(H,M, A, lmax, dopol = want_pol)
        !This automatically frees previous content of A, and returns new one
        
        call HealpixAlm2Power(A,P)
        call HealpixAlm_Free(A)
        !Note usually no need to free objects unless memory is short

        call HealpixPower_Write(P,lensed_cls_file)

        !Save lensed map
        call HealpixMap_Write(M, lensed_map_file, overwrite=.true.)
        
    end if

#ifdef MPIPIX
    call HealpixFree(H)
    call mpi_finalize(i)
#endif

#ifdef DEBUG
    write (*,*) 'End of program'
    pause
#endif
    end program SimLensCMB