// JavaScript Document
<div id="twitter_data"><span style="font-style:italic;font-size:12px;">Loading tweets</span></div>



var twitter_id = 'dj__simba';
var twitter_count = 5;
var twitter_elem_id = 'twitter_data';

$(document).ready(function() {
	var url = 'http://twitter.com/statuses/user_timeline/'
	+ twitter_id
	+ '.json?callback=twitterCallback&count='
	+ twitter_count;
	var script = document.createElement('script');
	script.setAttribute('src', url);
	document.body.appendChild(script);
});

function twitterCallback(obj){
  var html = "";
  for(var i=0 ; i<twitter_count && i != obj.length; i++){
	html += "<div class='item'>";
		html += "<div class='twitterIcon'><img src='" +  obj[i].user.profile_image_url + "'> </div>";
		html += "<div class='bubbleLeft'>";
			html += "<div class='subheading'><strong>Tweet</strong> from <strong>" + obj[i].user.name + "</strong></div>";
			html += "<div class='tweetMessage'>" + obj[i].text + "</div>";
		html += "</div>";
		html += "<div class='bubbleRight'>&nbsp;</div>";
	html += "</div>";
	html += "<br clear='all' />";
}
  document.getElementById(twitter_elem_id).innerHTML = html;
}
