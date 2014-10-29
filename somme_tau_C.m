 

    for k = 1:K
        offset = (k-1) * 3;
        x = offset + 1;
        y = offset + 2;
        z = offset + 3;
                
        if ( is_contact(sample, k) == 1 )
            tau_C_x = epsilon(k) * ( contact_coord(sample, y) * normal_vectors(sample, z) - ...
                contact_coord(sample, z) * normal_vectors(sample, y) )

            tau_C_y = epsilon(k) * ( contact_coord(sample, x) * normal_vectors(sample, z) - ...
                contact_coord(sample, z) * normal_vectors(sample, x) )
        end
    end
    