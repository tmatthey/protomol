#!/usr/bin/perl 

##################################################################
use strict;
use Tk;
use Tk::ROText;
use IO::Handle;
use Tk::Dialog;
use Tk::ErrorDialog;
use Tk::NoteBook;
use Tk::TextUndo;

my $protomol = "";
my $configuration = "";
my $outputFile = "";
my $options = "";
my $textOutput;
my $textEdit;
my $select;
my $textInfo ="";
my $exe;
$exe->{-out}    = \$textOutput;
$exe->{-finish} = 0;
$exe->{-command} = "true ";
my $selection;
my $margin = 15;
my $last = "";
my $miny = "";
my $maxy = "";
my $tminy = "";
my $tmaxy = "";
my $cross = "";
my $fromX = "";
my $toX = "";
my $lineY = "";
my $avg = "Avg=";
my $SD= "SD=";
my $modifier;
my $lastUpdate = time();
my $autoUpdate = "";
my @header = ();
my @data = ();
my $showAvg = 1;

if($^O eq 'MSWin32'){
    $modifier = 'Control';
    $protomol = "../Release/protomol.exe";
}
elsif ($^O eq 'MacOS') {  # one of these days
    $modifier = 'Command';
}
else {
    $protomol = "../applications/protomol-app/protomol";
    $modifier = 'Control';
}


sub kill_command {};
sub killfam {};
sub get_pids {};
sub execute_command {};
sub plotit {};

#
# Top Window
my $top = MainWindow->new();
$top->title("Perl/Tk ProtoMol");
$top->iconname("ProtoMol");
$top->geometry("+32+32");
$exe->{-top}    = \$top;


#
# Menu
my $menubar = $top->Menu(-type => 'menubar');
$top->configure(-menu => $menubar);
my $menuFile  = $menubar->cascade(-label => '~File', -tearoff => 0);
my $open;
my $reopen;
my $saveAs;
my $save;
my $quit; 

# New
$menuFile->command(-label => 'New', 
		   -underline => 0,
		   -command => sub  {
		       $textInfo = "";
		       my $popupFile = $textEdit->menu()->entrycget('File',-menu);
		       $popupFile->invoke($popupFile->index('Clear'))."\n";
		       $configuration = undef;
		       $textEdit->FileName($configuration);
		       $textInfo = "Text cleared."
		       });
$menuFile->separator;

# Open
$open = $menuFile->command(-label => 'Open ...',
			   -underline => 0,
			   -accelerator => "$modifier-o",
			   -command => sub {
			       $textInfo = "";
			       return if($textEdit->numberChanges() > 0 &&  $top->Dialog(-title => "Dialog",
											 -text => "The text has been modified without being saved.",
											 -default_button => "Cancel",
											 -buttons => ["Cancel", "Open"],
											 -bitmap => 'question')->Show() eq 'Cancel');
			       $textEdit->FileName(undef);
			       my $popupFile = $textEdit->menu()->entrycget('File',-menu);
			       $popupFile->invoke($popupFile->index('Open'))."\n";
			       $configuration = $textEdit->FileName();
			       $textInfo = "File \'".$configuration."\' opened." if(defined($configuration) && length($configuration));
			       changeDir($configuration);
			       plotit();
			   });
$top->bind("<$modifier-o>" => sub {my $c = $menubar->entrycget($menuFile->cget(-label), -menu);$c->invoke($c->index($open->cget(-label)));});

# Reopen
$reopen = $menuFile->command(-label => 'Reopen ...', 
			     -underline => 0,
			     -accelerator => "$modifier-r",
			     -command => sub  {
				 $textInfo = "";
				 return if($textEdit->numberChanges() > 0 &&  $top->Dialog(-title => "Dialog",
											   -text => "The text has been modified without being saved.",
											   -default_button => "Cancel",
											   -buttons => ["Cancel", "Open"],
											   -bitmap => 'question')->Show() eq 'Cancel');
				 my $popupFile = $textEdit->menu()->entrycget('File',-menu);
				 if($textEdit->FileName() eq undef){
				     $popupFile->invoke($popupFile->index('Open'))."\n";
				 }
				 else {
				     $textEdit->Load($configuration);
				 }
				 $configuration = $textEdit->FileName();
				 $textInfo = "File \'".$configuration."\' re-opened." if(defined($configuration) && length($configuration)) ;
				 changeDir($configuration);
				 plotit();
			     });
$top->bind("<$modifier-r>" => sub {my $c = $menubar->entrycget($menuFile->cget(-label), -menu);$c->invoke($c->index($reopen->cget(-label)));});

# Save
$menuFile->separator;
$save = $menuFile->command(-label => 'Save', 
			   -underline => 0,
			   -accelerator => "$modifier-s",
			   -command => sub  {
			       $textInfo = "";
			       my $popupFile = $textEdit->menu()->entrycget('File',-menu);
			       $popupFile->invoke($popupFile->index('Save'))."\n";
			       $configuration = $textEdit->FileName();
			       $textInfo = "File \'".$configuration."\' saved." if(defined($configuration) && length($configuration)) ;				 
			       changeDir($configuration);
			   });
$top->bind("<$modifier-s>" => sub {my $c = $menubar->entrycget($menuFile->cget(-label), -menu);$c->invoke($c->index($save->cget(-label)));});

# Save as
$saveAs = $menuFile->command(-label => 'Save As ...', 
			     -underline => 1,
			     -accelerator => "$modifier-a",
			     -command => sub  {
				 $textInfo = "";
				 my $popupFile = $textEdit->menu()->entrycget('File',-menu);
				 $popupFile->invoke($popupFile->index('Save As'))."\n";
				 $configuration = $textEdit->FileName();
				 $textInfo = "File \'".$configuration."\' saved as." if(defined($configuration) && length($configuration)) ;				 
				 changeDir($configuration);
			     });
$top->bind("<$modifier-a>" => sub {my $c = $menubar->entrycget($menuFile->cget(-label), -menu);$c->invoke($c->index($saveAs->cget(-label)));});

# Quit
$menuFile->separator;
$quit = $menuFile->command(-label => 'Quit', 
			   -underline => 0,
			   -accelerator => "$modifier+q",
			   -command => sub { 
			       if($textEdit->numberChanges() > 0){
				   my $answer = $top->Dialog(-title => "Quit Requester",
							     -text => "Actual file was modified or is unsaved.",
							     -default_button => "Save and Quit",
							     -buttons => ["Cancel", "Save As and Quit", "Save and Quit", "Quit"],
							     -bitmap => 'question')->Show();
				   
				   return if ($answer eq 'Cancel');
				   my $c = $menubar->entrycget($menuFile->cget(-label), -menu);
				   $c->invoke($c->index($save->cget(-label))) if ($answer eq "Save and Quit");
				   $c->invoke($c->index($saveAs->cget(-label)))    if ($answer eq "Save As and Quit");
			       }
			       kill_command;
			       exit 0;
			   });    
$top->bind("<$modifier-q>" => sub {my $c = $menubar->entrycget($menuFile->cget(-label), -menu);
				   $c->invoke($c->index($quit->cget(-label)));});


# Options
my $menuOptions  = $menubar->cascade(-label => '~Options', -tearoff => 0);

# Binary
$menuOptions->command(-label => 'ProtoMol binary ...', 
		      -command => sub {
			  my $fn = $top->getOpenFile(-parent   => $top,
						     -title    => 'ProtoMol Binary', 
						     -filetypes=>[['Binaries',       ['.exe', '.EXE', 'protomol'],],
								  ['All Files',        '*',             ],]);
			  $protomol = $fn if defined($fn) && length($fn);
			  $textInfo = "Binary \'".$protomol."\' opened.";
		      });

# Follow graph
$menuOptions->checkbutton(-label => "Follow graph", 
			  -variable => \$cross, -onvalue => "1",
			  -offvalue => "");

# Help  
my $menuHelp          = $menubar->cascade(-label => '~Help', -tearoff => 0);
$menuHelp->command(-label => 'About',
		   -command => sub {
		       $top->Dialog(-title => "TkProtoMol",
				    -text => "TkProtoMol",
				    -default_button => "Thanks",
				    -bitmap => 'info',
				    -buttons => ["Thanks"]
				    )->Show();
		   });

#
# textInfo
my $frame0 = $top->Frame()->pack(-expand => 'no', -fill => 'x',  -side => 'bottom');
$frame0->Label(-textvariable => \$textInfo)->pack(  -side => 'left');

#
# Notebook: Run
my $notebook = $top->NoteBook()->pack(-expand => 'yes', -fill => 'both',  -side => 'top');

# Run
my $nb1 = $notebook->add('Run', -label => 'Run');
my $frame1 = $nb1->Frame()->pack(-expand => 'no',-fill => 'x',  -side => 'top');

# Binary
$frame1->Label(-text => "ProtoMol")->grid( -row => 2, -column => 0, -sticky => 'ew');
$frame1->Entry( -textvariable  => \$protomol)->grid( -row => 2, -column => 1, -sticky => 'we');

# Configuration
$frame1->Label(-text => "Configuration")->grid( -row => 0, -column => 0, -sticky => 'ew');
$frame1->Entry(-textvariable  => \$configuration)->grid( -row => 0, -column => 1, -sticky => 'we');

# Options
$frame1->Label(-text => "Options")->grid( -row => 1, -column => 0, -sticky => 'ew');
$frame1->Entry(-textvariable  => \$options)->grid( -row => 1, -column => 1, -sticky => 'we');

$frame1->gridColumnconfigure(0, -weight =>0);
$frame1->gridColumnconfigure(1, -weight =>1);

# Run
my $frame2 = $nb1->Frame()->pack(-expand => 'no', -fill => 'x',  -side => 'top');
my $doit = $frame2->Button(-text => "Run",  
			   -command =>sub { 
			       if($textEdit->numberChanges() > 0 && 
				  $top->Dialog(-title => "Dialog",
					       -text => "The text has been modified without being saved.",
					       -default_button => "Save and Run",
					       -buttons => ["Save and Run", "Run"],
					       -bitmap => 'question')->Show() eq 'Save and Run'){
				   $textInfo = "";
				   my $popupFile = $textEdit->menu()->entrycget('File',-menu);
				   $popupFile->invoke($popupFile->index('Save'))."\n";
				   $configuration = $textEdit->FileName();
				   $textInfo = "File \'".$configuration."\' saved." if(defined($configuration) && length($configuration)) ;				 
				   changeDir($configuration);           
			       }
			       my $e1 = (length($protomol) > 0 ? "\"$protomol\"":$protomol);
			       my $e2 = (length($configuration) > 0 ? "\"$configuration\"":$configuration);
			       $exe->{-command} = "$e1 $e2 $options";
			       $exe->{-finish} = 0;
			       $textInfo = "Running \'".$configuration."\' ...";				 
			       execute_command;				     
			   })->grid( -row => 0, -column => 0, -sticky => 'ew');
$exe->{-doit} = \$doit;

# Clear
$frame2->Button(-text => "Clear",  
		-command =>sub { 
		    $textOutput->delete("1.0", "end");	
		})->grid( -row => 0, -column => 1, -sticky => 'ew');

$frame2->gridColumnconfigure(0, -weight =>1);
$frame2->gridColumnconfigure(1, -weight =>1);


# textOutput
$textOutput = $nb1->Scrolled('ROText', -scrollbars => 'se')->pack(-expand => 'yes', -fill => 'both',  -side => 'top');
$textOutput->configure(-wrap => 'none');
$textOutput->configure(-state => "normal");

#
# Notebook: Edit
my $nb2 = $notebook->add('Edit', -label => 'Edit');
my $frame3 = $nb2->Frame()->pack(-fill => 'x',  -side => 'top');

# Configuration
$frame3->Label(-text => "Configuration")->grid( -row => 0, -column => 0, -sticky => 'ew');
$frame3->Entry(-textvariable  => \$configuration)->grid( -row => 0, -column => 1, -sticky => 'we');
$frame3->gridColumnconfigure(0, -weight =>0);
$frame3->gridColumnconfigure(1, -weight =>1);

# textEdit
$textEdit = $nb2->Scrolled('TextUndo', -scrollbars => 'se')->pack(-expand => 'yes', -fill => 'both',  -side => 'top');
$textEdit->configure(-wrap => 'none');
my $popupFile = $textEdit->menu()->entrycget('File',-menu);
$textEdit->menu()->configure(-postcommand => sub{$textEdit->FileName($configuration)});

#
# Notebook : plot
my $nb3 = $notebook->add('Plot', -label => 'Plot');
my $frame4 = $nb3->Frame()->pack(-fill => 'x',  -side => 'top');

# plot
$frame4->Button(-text => "Plot",  -command => sub  {plotit();})->grid( -row => 0, -column => 0, -sticky => 'ew');

$frame4->Label(-text => "From:")->grid( -row => 0, -column => 1, -sticky => 'ew');
$frame4->Entry(-textvariable  => \$fromX,  -width => 6)->grid( -row => 0, -column => 2, -sticky => 'we');
$frame4->Label(-text => "To:")->grid( -row => 0, -column => 3, -sticky => 'ew');
$frame4->Entry(-textvariable  => \$toX,  -width => 6)->grid( -row => 0, -column => 4, -sticky => 'we');
$frame4->Label(-text => "Line:")->grid( -row => 0, -column => 5, -sticky => 'ew');
$frame4->Entry(-textvariable  => \$lineY,  -width => 6)->grid( -row => 0, -column => 6, -sticky => 'we');

# Autoupdate
$frame4->Checkbutton(-text => "Auto update", 
		     -variable => \$autoUpdate, 
		     -onvalue => "1",
		     -offvalue => "")->grid( -row => 0, -column => 7, -sticky => 'we');

$frame4->Entry(-textvariable  => \$outputFile)->grid( -row => 0, -column => 8, -sticky => 'we');

# browse
$frame4->Button(-text => "Browse ...",
		-command => sub {
		    my $fn = $top->getOpenFile(-parent   => $top,
					       -title    => 'Output File', 
					       -filetypes=>[['All Files',        '*',             ],]);
		    $outputFile = $fn if defined($fn) && length($fn);
		    $textInfo = "Outpufile \'".$outputFile."\' opened.";
		    plotit();
		})->grid( -row => 0, -column => 9, -sticky => 'ew');
$frame4->gridColumnconfigure(8, -weight =>1);
my $frame5 = $nb3->Frame()->pack(-expand => 'yes', -fill => 'both',  -side => 'top');

# canvas plot
my $plot = $frame5->Canvas(-relief => 'groove',
			   -borderwidth => '5'
			   )->pack(-side => 'right', -expand => 'yes', -fill => 'both');
$plot->Tk::bind('<Motion>' => \&display_coordinates);
$plot->Tk::bind('<Button-1>' => sub{
    my($canvas) = @_;

    my $e = $canvas->XEvent;
    my($x, $y) = ($e->x, $e->y); 
	

    my $minx = ($fromX eq ""?0:$fromX);
    my $maxx = ($toX eq ""?$#data:$toX);
    $minx = 0      if($minx < 0);
    $maxx = $#data if($#data < $maxx);
    if($minx>=$maxx){
	$minx = 0;
	$maxx = $#data;
    }
    my $w = $plot->width();
    my $h = $plot->height();
    if($x <= $margin){
	$fromX = "";
    }
    elsif($x >= $w-$margin){
	$fromX = "";
	$toX = "";
    }
    else{
	$fromX=$minx+int(($x-$margin)*($maxx+1-$minx)/($w-2.0*$margin));
    }
    plotit();

});
$plot->Tk::bind('<Button-3>' => sub{
    my($canvas) = @_;

    my $e = $canvas->XEvent;
    my($x, $y) = ($e->x, $e->y); 
	

    my $minx = ($fromX eq ""?0:$fromX);
    my $maxx = ($toX eq ""?$#data:$toX);
    $minx = 0      if($minx < 0);
    $maxx = $#data if($#data < $maxx);
    if($minx>=$maxx){
	$minx = 0;
	$maxx = $#data;
    }
    my $w = $plot->width();
    my $h = $plot->height();
    if($x <= $margin){
	$toX = "";
	$fromX = "";
    }
    elsif($x >= $w-$margin){
	$toX = "";
    }
    else{
	$toX=$minx+int(($x-$margin)*($maxx+1-$minx)/($w-2.0*$margin));
    }
    plotit();

});
#$plot->Tk::bind('<Configure>' => sub { plotit() if($outputFile && $lastUpdate< time());});



# bind quit
$top->protocol('WM_DELETE_WINDOW' => sub { my $c = $menubar->entrycget($menuFile->cget(-label), -menu);
					   $c->invoke($c->index($quit->cget(-label)));});



# bind wheel
if ($^O eq 'MSWin32') {
    $top->bind('<MouseWheel>' =>
	       [ sub { $textEdit->yview('scroll', -($_[1] / 120) * 3, 'units') },
		 Ev('D') ]
	       );
}
else {
    
    # Support for mousewheels on Linux commonly comes through
    # mapping the wheel to buttons 4 and 5.  If you have a
    # mousewheel ensure that the mouse protocol is set to
    # "IMPS/2" in your /etc/X11/XF86Config (or XF86Config-4)
    # file:
    #
    # Section "InputDevice"
    #     Identifier  "Mouse0"
    #     Driver      "mouse"
    #     Option      "Device" "/dev/mouse"
    #     Option      "Protocol" "IMPS/2"
    #     Option      "Emulate3Buttons" "off"
    #     Option      "ZAxisMapping" "4 5"
    # EndSection	
    $top->bind('<4>' => sub {
	$textEdit->yview('scroll', -3, 'units') unless $Tk::strictMotif;
    });
    
    $top->bind('<5>' => sub {
	$textEdit->yview('scroll', +3, 'units') unless $Tk::strictMotif;
    });
}

MainLoop;

exit;

##################################################################
sub display_coordinates {
##################################################################

    # Print Canvas and Plot coordinates.

    my($canvas) = @_;

    my $e = $canvas->XEvent;
    my($x, $y) = ($e->x, $e->y); 
    $plot->delete('cross');
    if($#data < 0 || $data[0] =~ /^\s*$/){
	plotit();
    }
    if($#data < 0 || $data[0] =~ /^\s*$/){
	$textInfo = "Empty file \'".$outputFile."\'.";
	return;
    }
	

    my $minx = ($fromX eq ""?0:$fromX);
    my $maxx = ($toX eq ""?$#data:$toX);
    $minx = 0      if($minx < 0);
    $maxx = $#data if($#data < $maxx);
    if($minx>=$maxx){
	$minx = 0;
	$maxx = $#data;
    }
    my $w = $plot->width();
    my $h = $plot->height();
    my $i=0;
    if($x <= $margin){
	$i = $minx;
    }
    elsif($x >= $w-$margin){
	$i = $maxx;
    }
    else{
	$i=$minx+int(($x-$margin)*($maxx+1-$minx)/($w-2.0*$margin));
    }
    my @a0 = split(/\s+/,$data[$i]);
    $plot->create('line',$x,$margin,$x,$h-$margin,
		  -tags => ['cross'], -fill => 'blue') if($x>=$margin && $x <= $w-$margin);
    my $y2 = "";
    if($cross){
	my $ay = ($h-2*$margin)/($miny < $maxy ? $maxy-$miny: 1.0);
	my $y1 = $a0[$selection]-$miny;  
	$y1 = $h-($margin+$y1*$ay);
	$plot->create('line',$margin,$y1,$w-$margin,$y1,
		      -tags => ['cross'], -fill => 'blue') if($y1>=$margin && $y1 <= $h-$margin);
    }
    else {
	if($y>=$margin && $y <= $h-$margin){
	    $plot->create('line',$margin,$y,$w-$margin,$y,
			  -tags => ['cross'], -fill => 'blue') ;
	    $y2 =  " (".(($h-$y-$margin)*($maxy-$miny)/($h-2*$margin)+$miny).")";
	}
	
    }
    $textInfo = $header[$selection].": ".$a0[$selection].$y2." (i=$i)";
}

##################################################################
sub plotit {
##################################################################
    
    if(!$outputFile && $configuration){
	if(open(FILE,"<$configuration")){
	    my @config=<FILE>;
	    close(FILE);
	    chomp(@config);
	    # grab last outputfile
	    foreach my $i (0 .. $#config){		
		$config[$i] =~ s/^\s+//; 
		$config[$i] =~ s/\#.*$//; 
		my @tmp = split(/\s+/,$config[$i]);
		if(lc($tmp[0]) =~ /file$/ && !(lc($tmp[0]) =~ /^fin/) && $tmp[1] 
		   && lc($tmp[0]) ne 'posfile' 
		   && lc($tmp[0]) ne 'parfile' 
		   && lc($tmp[0]) ne 'psffile' ){
		    my $fn = lc($tmp[1]);
		    $outputFile = $tmp[1] if($fn ne 'no' && $fn ne 'yes' && $fn ne 'true' && $fn ne 'false');
		}
	    }
	    if($outputFile){
		my $str = $configuration;
		$str =~ s/[^\\\/]*$//g;
		$str =~ s/^\s+//g;
		$str =~ s/\s+$//g;
		$outputFile = $str.$outputFile;
	    }    
	}
    }
    @data = ();
    @header = ();
    $plot->delete('plot');
    $plot->delete('cross');
    if(open(FILE,"<$outputFile")){
	@data=<FILE>;
	close(FILE);
	chomp(@data);
	foreach my $i (0 .. $#data){
	    $data[$i] =~ s/^\s+//; 
	}
	if($#data < 0 || $data[0] =~ /^\s*$/){
	    $textInfo = "Empty file \'".$outputFile."\'.";
	    return;
	}	
	if(open(FILEH,"<$outputFile.header")){
	    @header=<FILEH>;
	    chomp(@header);
	    $header[0] =~ s/^\s+//;
	    @header = split(/\s+/,$header[0]);
	    close(FILEH);	    
	}
	else {
	    @header = split(/\s+/,$data[0]);
	    foreach my $i (0 .. $#header){
		$header[$i] = $i+1;
	    }
	    $textInfo = "Could not find header of \'".$outputFile."\'.";	    
	}
    }
    else {
	$textInfo = "Can't open \'".$outputFile."\'.";
	return;
    }
    if(!$select || $last ne $outputFile){
	$select->destroy() if($select);
	$select = $frame5->Frame()->pack(-side => 'left', -fill => 'y');
	$select->Label(-textvariable => \$tmaxy)->grid( -row => 0, -column => 0, -sticky => 'n');
	foreach my $i (0 .. $#header){
	    $select->Radiobutton( -text  => "$header[$i]",
				  -variable => \$selection,
				  -justify => 'left',
				  -value => $i,
				  -command => \&plotit
				  )->grid( -row => $i+1, -column => 0, -sticky => 'wn');
	}
	$select->Checkbutton(-variable => \$showAvg, 
			     -onvalue => "1",
			     -offvalue => "",
			     -command => \&plotit,
			     -textvariable => \$avg)->grid( -row => $#header+2, -column => 0, -sticky => 'n');
	$select->Label(-textvariable => \$SD)->grid( -row => $#header+3, -column => 0, -sticky => 'n');
	$select->Label(-textvariable => \$tminy)->grid( -row => $#header+4, -column => 0, -sticky => 's');
	$select->gridRowconfigure(0, -weight =>1);
	$select->gridRowconfigure($#header+4, -weight =>1);
	$selection = 0;
	$last = $outputFile;
    }
#    $top->update();
    my $minx = ($fromX eq ""?0:$fromX);
    my $maxx = ($toX eq ""?$#data:$toX);
    $minx = 0      if($minx < 0);
    $maxx = $#data if($#data < $maxx);
    if($minx>=$maxx){
	$minx = 0;
	$maxx = $#data;
    }
    foreach my $i ($minx .. $maxx){
	my @a0 = split(/\s+/,$data[$i]);
	$miny = $a0[$selection] if($miny > $a0[$selection] || $i == $minx);
	$maxy = $a0[$selection] if($maxy < $a0[$selection] || $i == $minx);
    }
    $tmaxy = "Max=".sprintf("%f",$maxy);
    $tminy = "Min=".sprintf("%f",$miny);
    my $w = $plot->width();
    my $h = $plot->height();
    my $ax = ($w-2*$margin)/($minx < $maxx ? $maxx-$minx: 1.0);
    my $ay = ($h-2*$margin)/($miny < $maxy ? $maxy-$miny: 1.0);
    my $sum = 0.0;
    my $n = 0;
    my @a1 =(); 
    my @a2 =();
    my $x1 =""; 
    my $y1 =""; 
    my $x2 =""; 
    my $y2 =""; 
    foreach my $i ($minx .. $maxx){
	@a1 = split(/\s+/,$data[$i]);
	$x1 = $i-$minx;
	$y1 = $a1[$selection]-$miny;  
	$x1 = $margin+$x1*$ax;
	$y1 = $margin+$y1*$ay;
	$sum += $a1[$selection];	
	$n++;
	$plot->create('line',$x1,$h-$y1,$x2,$h-$y2,
		      -tags => ['plot']) unless($x2 eq "");
	@a2 = @a1;
	$x2 = $x1;
	$y2 = $y1;
    }   
    if($lineY ne ""){
	$y1 = $lineY-$miny;  
	$y1 = $h-($margin+$y1*$ay);
	$plot->create('line',$margin,$y1,$w-$margin,$y1,
		     -tags => ['plot'], -fill => 'red');
    }
    if($n > 0){
	$avg = "Avg=".sprintf("%f",$sum/$n);
	if($showAvg){
	    $y1 = $sum/$n-$miny;  
	    $y1 = $h-($margin+$y1*$ay);
	    $x1 = $minx-$minx;
	    $x1 = $margin+$x1*$ax;
	    $x2 = $maxx-$minx;
	    $x2 = $margin+$x2*$ax;
	    $plot->create('line',$x1,$y1,$x2,$y1,
			  -tags => ['plot'], -fill => 'green');
	    $plot->create('line',$x1,$h-$margin,$x1,$margin,
			  -tags => ['plot'], -fill => 'green');
	    $plot->create('line',$x2,$h-$margin,$x2,$margin,
			  -tags => ['plot'], -fill => 'green');
	}
	my $tmp = $sum/$n;
	$sum = 0.0;
	if($n > 1){
	    foreach my $i ($minx .. $maxx){
		my @a1 = split(/\s+/,$data[$i]);
		my $x1 = $a1[$selection];  
		$sum += ($x1-$tmp)*($x1-$tmp);
	    }
	    $SD = "SD=".sprintf("%f",sqrt($sum/($n-1)));
	}
    }
    $textInfo = "plot(0,".($#data).").";
    $lastUpdate = time();
}


##################################################################
sub changeDir {
##################################################################
    my $str = shift;
    $str =~ s/[^\\\/]*$//g;
    $str =~ s/^\s+//g;
    $str =~ s/\s+$//g;
    chdir "\"".$str."\"" if(length($str)>0);
}

##################################################################
sub removeDirPath {
##################################################################
    my $str = shift;
    $str =~ s/^.*[\\\/]//g;
    $str =~ s/^\s+//g;
    $str =~ s/\s+$//g;
    return $str;
}

##################################################################
sub read_stdout {
##################################################################

    # Called when input is available for the output window.  Also checks
    # to see if the user has clicked Cancel.

    if ($exe->{-finish}) {
        kill_command;
    } else {
        my $h = $exe->{-handle};
        if ( sysread $h, $_, 4096 ) {
	    my $out = $_;
	    $out =~ s/\r+//g;	# Remove possible \r for Windows    
            my $t = ${$exe->{-out}};
            $t->insert('end', $out);
            $t->yview('end');
	    $nb1->idletasks();
        } else {
            $exe->{-finish} = 1;
        }
	plotit() if($autoUpdate && $lastUpdate+2 < time());
    }	
}

##################################################################
sub reset_doit_button {
##################################################################

    # Establish normal "Do It" button parameters.

    ${$exe->{-doit}}->configure(-text       => 'Run',
				-relief     => 'raised',
				-state      => 'normal',
				-command    => sub {$exe->{-finish} = 0;
						    if($textEdit->numberChanges() > 0 && 
						       $top->Dialog(-title => "Dialog",
								    -text => "The text has been modified without being saved.",
								    -default_button => "Save and Run",
								    -buttons => ["Save and Run", "Run"],
								    -bitmap => 'question')->Show() eq 'Save and Run')
						    {
							$textInfo = "";
							my $popupFile = $textEdit->menu()->entrycget('File',-menu);
							$popupFile->invoke($popupFile->index('Save'))."\n";
							$configuration = $textEdit->FileName();
							$textInfo = "File \'".$configuration."\' saved." if(defined($configuration) && length($configuration)) ;				 
							changeDir($configuration);           
						    }
						    my $e1 = (length($protomol) > 0 ? "\"$protomol\"":$protomol);
						    my $e2 = (length($configuration) > 0 ? "\"$configuration\"":$configuration);
						    $exe->{-command} = "$e1 $e2 $options";
						    $exe->{-finish} = 0;
						    execute_command;				     
						}
						
				);
   		       $textInfo = "Done."
}

##################################################################
sub execute_command {
##################################################################

    # Execute the command and capture stdout/stderr.

    my $h = IO::Handle->new;
    die "IO::Handle->new failed." unless defined $h;
    $exe->{-handle} = $h;
    ${$exe->{-out}}->delete("1.0", "end");	
    $exe->{-pid} = open $h, $exe->{-command} . ' 2>&1 |';
    if (not defined $exe->{-pid}) {
        ${$exe->{-out}}->insert('end', "Sorry, '" . $exe->{-command} . "' : $!\n");
        kill_command;
        return;
    }
    $h->autoflush(1);
    ${$exe->{-top}}->fileevent($h, 'readable' => [\&read_stdout]);

    ${$exe->{-doit}}->configure(
        -text    => 'Cancel',
        -relief  => 'raised',
        -state   => 'normal',
        -command => [\&kill_command],
    );

    
}


##################################################################
sub kill_command {
##################################################################
    
    # A click on the blinking Cancel button resumes normal operations.

    $exe->{-finish} = 1;
    my $h = $exe->{-handle};
    return unless defined $h;
    ${$exe->{-top}}->fileevent($h, 'readable' => ''); # clear handler
    killfam 'TERM', $exe->{-pid} if defined $exe->{-pid};
    close $h;
    reset_doit_button;
}

##################################################################
sub killfam {
##################################################################

    my($signal, @pids) = @_;

    eval "require Proc::ProcessTable";
    if (!$@) {
	my $pt = Proc::ProcessTable->new;
	my(@procs) =  @{$pt->table};
	my(@kids) = get_pids \@procs, @pids;
	@pids = (@pids, @kids);
    }

    kill $signal, @pids;

} # end killfam

##################################################################
sub get_pids {
##################################################################

    my($procs, @kids) = @_;

    my @pids;
    foreach my $kid (@kids) {
	foreach my $proc (@$procs) {
	    if ($proc->ppid == $kid) {
		my $pid = $proc->pid;
		push @pids, $pid, get_pids $procs, $pid;
	    } 
	}
    }
    @pids;

} # end get_pids
