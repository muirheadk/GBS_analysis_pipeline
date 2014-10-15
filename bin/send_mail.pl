=head1 NAME
 
 B<extract_ipod_music.pl> - Extract music from your iPod.
 
 =head1 DESCRIPTION
 
 This program copies all the music from your iPod to your desktop
 or directory of choice using the MP4::Info Perl module found on CPAN.
 Download using the following command 'sudo ./install.sh -v' as the program
 is dependant on several modules. (-v is for VERBOSE)
 
=cut

#!/usr/bin/perl
use warnings;
use strict;

use Mail::Mailer;
#use Switch;
#use POSIX;

#use constant DATE => strftime("%Y-%m-%d", localtime);
#use constant TIME => strftime("%%H-%M-%S", localtime);

use Getopt::Long;
my ($subject, $message, $mailer_type);
my @options = (
    's=s'    => \$subject,
    'm=s'    => \$message,
    't=s'    => \$mailer_type,
);
&GetOptions(@options);

usage() unless (
    defined $subject
    and defined $message
);
$mailer_type = 'email' unless defined $mailer_type;

my $to_address = 'kevin5@ualberta.ca';
my $from_address = 'cookeadmin <root@cookelab.ccis.ualberta.ca>';
my $email_to_text_address = '14039195832@pcs.rogers.com';

sub usage {
    
    die <<"USAGE";
    
Usage: $0 -s subject -m message -t mailer_type
    
    
    -s subject - 
    -m message - 
    -t mailer_type -  

USAGE
}



# Send email to the following address.
send_mail($to_address,$from_address,$subject,$message) if($mailer_type eq 'email' or $mailer_type eq 'both');

# Send email using email to text from your cellular provider.
send_mail($email_to_text_address,$from_address,$subject,$message) if($mailer_type eq 'text' or $mailer_type eq 'both');

=head1 SUBROUTINES
 
 $output_dir = get_output_dir($input_dir) - Generate iPod output directory using the specified iPod input directory.
 
 Input paramater(s):
 
 $input_dir - The iPod input directory.
 
 Output paramater(s):
 
 $output_dir - The iPod output directory.
 
=cut
sub send_mail{
    my $to_address = shift or die "lost recipient email address";
    my $from_address = shift or die "lost sender email address";
    my $subject = shift or die "lost email subject";
    my $message = shift or die "lost email message";
    
    warn "Sending mail with $subject to $to_address\n";
    my $mailer = new Mail::Mailer("sendmail");
    $mailer->open( {
        To       => $to_address,
        From     => $from_address,
        Subject  => $subject
    } );
    
    print $mailer <<END_OF_MESSAGE;
    $message
END_OF_MESSAGE
        
    close $mailer;
}

