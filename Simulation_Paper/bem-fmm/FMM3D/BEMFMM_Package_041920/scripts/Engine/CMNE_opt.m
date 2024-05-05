function MNE_out = CMNE_opt(MNE_in,P_tmp,t_tmp)
%%% define the maximum allocated current to each coil in this script
Gx = MNE_in.Gxu; 
Gy = MNE_in.Gyu; 
Gz = MNE_in.Gzu; 

            G = cat(1,Gx,Gy,Gz);

            [U S V]=svd(G,'econ');
            n_cutoff = 14;
            S_tmp_diag = diag(S);
            
            anal_all{4,1,1}.G_cond = S_tmp_diag(1)/S_tmp_diag(end);
            
            sc_tmp=1;

            lambda=sc_tmp*S_tmp_diag(n_cutoff);
         
            G_pinv=(G'*G+lambda^2*eye(size(G,2)))\G';
            %disp('Staring MNE...');
            

            dI_dt = 1e6; 
            P_tmp_x=zeros(length(Gx),1);
            P_tmp_y=zeros(length(Gy),1);
            P_tmp_z=zeros(length(Gz),1);
            
            Lia = ismember(t_tmp,MNE_in.targ_loc);
            Lia2 = sum(Lia,2);
            Ed_ind = find(Lia2 >=1);
            
            P_tmp_x(Ed_ind(1))=squeeze(MNE_in.targ_dirs(1,1,1));
            P_tmp_y(Ed_ind(2))=squeeze(MNE_in.targ_dirs(1,1,2));
            P_tmp_z(Ed_ind(3))=squeeze(MNE_in.targ_dirs(1,1,3));
 
            dI_dt_MNE=G_pinv*(cat(1,P_tmp_x,P_tmp_y,P_tmp_z));     
            scale_dI_dt=dI_dt/max(abs(dI_dt_MNE));
            dI_dt_MNE=scale_dI_dt*dI_dt_MNE;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%% CONSTRAINED MNE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
             options = optimset('DerivativeCheck','off','GradObj','on','Display','off');
             E_max_targ=100;
             nq=16;
             n_coils=16;
             max_current=80*10^(6);
        
            %dI_dt_MNE_con_weig=zeros(n_targ_dirs,length(inds_array_th),n_coils);

            %%%%%%%%%%%% THE SVD MAKES IT CONVENIENT TO MODEL
            %%%%%%%%%%%% REALIZABLE E-FIELDS   
        
            E_x_tmp=Gx*dI_dt_MNE;
            E_y_tmp=Gy*dI_dt_MNE;
            E_z_tmp=Gz*dI_dt_MNE;
            E_abs_tmp=sqrt(E_x_tmp.^2+E_y_tmp.^2+E_z_tmp.^2);
    
            dI_dt_MNE_tmp= dI_dt_MNE;
            
            %%%%SETTING THE MNE AS THE TARGET FIELD
            
            P_targ=cat(1,E_x_tmp,E_y_tmp,E_z_tmp); 
            unshimmed_vec=(P_targ/max(E_abs_tmp))*E_max_targ;
             
            %init=dI_dt_vals/max(abs(dI_dt_vals))*max_current;
     
            %%%INITIALIZE WITH SCALED SOLUTIONS
            init=dI_dt_MNE_tmp/max(abs(dI_dt_MNE_tmp))*max_current;
        
            %amps = fmincon( @(amps)shim_fun(amps,A,unshimmed_vec),zeros(nq,1),[],[],[],[],-max_current*ones(nq,1),max_current*ones(nq,1),[],options);
            amps = fmincon( @(amps)shim_fun(amps,G,unshimmed_vec),init,[],[],[],[],-max_current*ones(nq,1),max_current*ones(nq,1),[],options);
             
            dI_dt_MNE_con=amps;
            
            E_x_con=Gx*dI_dt_MNE_con;
            E_y_con=Gy*dI_dt_MNE_con;
            E_z_con=Gz*dI_dt_MNE_con;
            E_abs_con=sqrt(E_x_con.^2+E_y_con.^2+E_z_con.^2);


    H_2=patch('Faces',t_tmp,'Vertices',P_tmp,'FaceVertexCData',E_abs_con,'FaceColor','flat','Edgecolor','none');hold on; 
    axis on; lighting gouraud; daspect([ 1 1 1]); 
    set(gcf,'Color','White'); colormap jet; h1 = colorbar; set(get(h1,'title'),'string','E (V/m)');
    view(-180,90);
    
    MNE_out.E_abs_con = E_abs_con;
    MNE_out.dI_dt = dI_dt_MNE_con;
end