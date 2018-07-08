	merge.with.order <- function(x,y, ..., sort = T, keep_order)
	{
	  
		# this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
		add.id.column.to.data <- function(DATA)
		{
			data.frame(DATA, id... = seq_len(nrow(DATA)))
		}
		# add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
		order.by.id...and.remove.it <- function(DATA)
		{
			# gets in a data.frame with the "id..." column.  Orders by it and returns it
			if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")
 
			ss_r <- order(DATA$id...)
			ss_c <- colnames(DATA) != "id..."
			DATA[ss_r, ss_c]		
		}
 
		# tmp <- function(x) x==1; 1	# why we must check what to do if it is missing or not...
		# tmp()
 
		if(!missing(keep_order))
		{
			if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
			if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
			# if you didn't get "return" by now - issue a warning.
			warning("The function merge.with.order only accepts NULL/1/2 values for the keep_order variable")
		} else {return(merge(x=x,y=y,..., sort = sort))}
	}
 
