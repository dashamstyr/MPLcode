import numpy as np

def hist_timealt(times,altitudes,vals,bin_times,bin_alts):

    time_count,time_index,lowtimes,hightimes=bin_times.do_bins(times)    
    alt_count,alt_index,lowalts,highalts=bin_alts.do_bins(alts)    
    out_vals=np.empty([bin_times.numbins,bin_alts.numbins],dtype=np.object)
    for row in range(bin_times.numbins):
        for col in range(bin_alts.numbins):
            out_vals[row,col]=list()

    for data_index in range(time_index.size):
        grid_row=time_index[data_index]
        grid_col=alt_index[data_index]
        if grid_row < 0 or grid_col < 0:
            continue
        try:
            out_vals[grid_row,grid_col].append(data_index)
        except:
            print "trouble at: data_index {0:d}, gives {1:d},{2:d}".format(data_index,grid_row,grid_col)

    val_grid=np.empty_like(out_vals,dtype=np.float)
    time_grid=np.empty_like(out_vals,dtype=np.float)
    alt_grid=np.empty_like(out_vals,dtype=np.float)
    rows,cols=val_grid.shape
    vals=vals.ravel()
    for the_row in range(rows):
        for the_col in range(cols):
            val_list=out_vals[the_row,the_col]
            if len(val_list)==0:
                val_grid[the_row,the_col]=np.nan
                time_grid[the_row,the_col]=bin_times.get_centers()[the_row]
                alt_grid[the_row,the_col]=bin_alts.get_centers()[the_col]
            else:
                try:
                    val_vals=np.take(vals,val_list)
                    time_vals=np.take(times,val_list)
                    alt_vals=np.take(alts,val_list)
                    val_grid[the_row,the_col]=np.sum(val_vals)
                    time_grid[the_row,the_col]=np.mean(time_vals)                    
                    alt_grid[the_row,the_col]=np.mean(alt_vals)                    
                except IndexError:
                    print "oops: ",val_list
    return time_grid,alt_grid,val_grid


class hist_2D(object):

    def __init__(self,times,alts,vals,bin_times,bin_alts):
        self.times=times
        self.alts=alts
        self.vals=vals
        self.bin_times=bin_times
        self.bin_alts=bin_alts
        self.time_count,self.time_index,self.lowtimes,self.hightimes=bin_times.do_bins(times)    
        self.alt_count,self.alt_index,self.lowalts,self.highalts=bin_alts.do_bins(alts)    
        self.out_vals=np.empty([bin_times.numbins,bin_alts.numbins],dtype=np.object)

    def calc_vals(self):

        numtimebins=self.bin_times.numbins
        numaltbins=self.bin_alts.numbins

        for row in range(numtimebins):
            for col in range(numaltbins):
                self.out_vals[row,col]=list()
        num_datapts=self.times.size
        for data_index in range(num_datapts):
            grid_row=self.time_index[data_index]
            grid_col=self.alt_index[data_index]
            if grid_row < 0 or grid_col < 0:
                continue
            else:
                self.out_vals[grid_row,grid_col].append(data_index)

    def calc_mean(self):
        
        self.calc_vals()
        val_grid=np.empty_like(self.out_vals,dtype=np.float32)
        time_grid=np.empty_like(self.out_vals,dtype=np.float32)
        alt_grid=np.empty_like(self.out_vals,dtype=np.float32)
        rows,cols=self.out_vals.shape
        
        for row in range(rows):
            for col in range(cols):
                val_list=self.out_vals[row,col]
                if len(val_list)==0:
                    val_grid[row,col]=np.nan
                    time_grid[row,col]=self.bin_times.get_centers()[row]
                    alt_grid[row,col]=self.bin_alts.get_centers()[col]
                else:
                    val_vals=np.take(self.vals,val_list)
                    time_vals=np.take(self.times,val_list)
                    alt_vals=np.take(self.alts,val_list)
                    val_grid[row,col]=np.mean(val_vals)
                    time_grid[row,col]=np.mean(time_vals)                    
                    alt_grid[row,col]=np.mean(alt_vals)                    

        return np.asarray(time_grid),np.asarray(alt_grid),np.asarray(val_grid)

    def calc_sum(self):
        
        self.calc_vals()
        val_grid=np.empty_like(self.out_vals,dtype=np.float32)
        time_grid=np.empty_like(self.out_vals,dtype=np.float32)
        alt_grid=np.empty_like(self.out_vals,dtype=np.float32)
        rows,cols=self.out_vals.shape
        
        for row in range(rows):
            for col in range(cols):
                val_list=self.out_vals[row,col]
                if len(val_list)==0:
                    val_grid[row,col]=np.nan
                    time_grid[row,col]=self.bin_times.get_centers()[row]
                    alt_grid[row,col]=self.bin_alts.get_centers()[col]
                else:
                    val_vals=np.take(self.vals,val_list)
                    time_vals=np.take(self.times,val_list)
                    alt_vals=np.take(self.alts,val_list)
                    val_grid[row,col]=np.sum(val_vals)
                    time_grid[row,col]=np.mean(time_vals)                    
                    alt_grid[row,col]=np.mean(alt_vals)                    

        return np.asarray(time_grid),np.asarray(alt_grid),np.asarray(val_grid)


class hist_class(object):

    def __init__(self,times,alts,vals,bin_times,bin_alts):
        self.times=times
        self.alts=alts
        self.vals=vals
        self.bin_times=bin_times
        self.bin_alts=bin_alts
        self.time_count,self.time_index,self.lowtimes,self.hightimes=bin_times.do_bins(times)    
        self.alt_count,self.alt_index,self.lowalts,self.highalts=bin_alts.do_bins(alts)    
        self.out_vals=np.empty([bin_times.numbins,bin_alts.numbins],dtype=np.object)

    def calc_vals(self):

        numtimebins=self.bin_times.numbins
        numaltbins=self.bin_alts.numbins

        for row in range(numtimebins):
            for col in range(numaltbins):
                self.out_vals[row,col]=list()
        num_datapts=self.times.size
        for data_index in range(num_datapts):
            grid_row=self.time_index[data_index]
            grid_col=self.alt_index[data_index]
            if grid_row < 0 or grid_col < 0:
                continue
            else:
                self.out_vals[grid_row,grid_col].append(data_index)

    def calc_mean(self):
        
        self.calc_vals()
        val_grid=np.empty_like(self.out_vals,dtype=np.float32)
        time_grid=np.empty_like(self.out_vals,dtype=np.float32)
        alt_grid=np.empty_like(self.out_vals,dtype=np.float32)
        rows,cols=self.out_vals.shape
        
        for row in range(rows):
            for col in range(cols):
                val_list=self.out_vals[row,col]
                if len(val_list)==0:
                    val_grid[row,col]=np.nan
                    time_grid[row,col]=self.bin_times.get_centers()[row]
                    alt_grid[row,col]=self.bin_alts.get_centers()[col]
                else:
                    val_vals=np.take(self.vals,val_list)
                    time_vals=np.take(self.times,val_list)
                    alt_vals=np.take(self.alts,val_list)
                    val_grid[row,col]=np.mean(val_vals)
                    time_grid[row,col]=np.mean(time_vals)                    
                    alt_grid[row,col]=np.mean(alt_vals)                    

        return np.asarray(time_grid),np.asarray(alt_grid),np.asarray(val_grid)

    def calc_sum(self):
        
        self.calc_vals()
        val_grid=np.empty_like(self.out_vals,dtype=np.float32)
        time_grid=np.empty_like(self.out_vals,dtype=np.float32)
        alt_grid=np.empty_like(self.out_vals,dtype=np.float32)
        rows,cols=self.out_vals.shape
        
        for row in range(rows):
            for col in range(cols):
                val_list=self.out_vals[row,col]
                if len(val_list)==0:
                    val_grid[row,col]=np.nan
                    time_grid[row,col]=self.bin_times.get_centers()[row]
                    alt_grid[row,col]=self.bin_alts.get_centers()[col]
                else:
                    val_vals=np.take(self.vals,val_list)
                    time_vals=np.take(self.times,val_list)
                    alt_vals=np.take(self.alts,val_list)
                    val_grid[row,col]=np.sum(val_vals)
                    time_grid[row,col]=np.mean(time_vals)                    
                    alt_grid[row,col]=np.mean(alt_vals)                    

        return np.asarray(time_grid),np.asarray(alt_grid),np.asarray(val_grid)