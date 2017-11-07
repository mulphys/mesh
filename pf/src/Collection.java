public class Collection<Item> {
	class List {
		Item item;
		List next,prev;
		List() {
			item=null;
			next=prev=null;
		}
		public void set(Item item) {
			this.item=item;
		}
	}
	final int ntags=2;
	int n;
	List current,first,tag,tag1;// tag[] did not work: "generic array creation" error!
	Collection() {
		n=0;
		first=current=tag=tag1=null;
	//?	tag=new List[ntags]; /// did not work: "generic array creation" error
	}
	public int number() { return n; }
	public int items() { return n; }
	public void add(Item item) {
		if(n==0){
			first=new List();
			first.item=item;
			first.prev=first.next=first;
			current=first;
		} else {//Add before the first:
			List newelement=new List();
			newelement.item=item;
			newelement.next=first;
			newelement.prev=first.prev;
			first.prev.next=newelement;
			first.prev=newelement;
		}
		n++;
	}
	public void insert(Item item) {
		if(n==0){
			first=new List();
			first.item=item;
			first.prev=first.next=first;
			current=first;
		} else {//Add after the current:
			List newelement=new List();
			newelement.item=item;
			newelement.next=current;
			newelement.prev=current.prev;
			current.prev.next=newelement;
			current.prev=newelement;
		}
		n++;
	}
	public void delete() {
		if(n==0) return;
		if(current==null) {
			System.out.println("CAN'T DELETE VERTEX");
			System.exit(1);
		}
		current.item=null;
		if(current==tag)tag=null;
		if(current==first)first=first.next;
		current.prev.next=current.next;
		current.next.prev=current.prev;
		current=current.next;
		if(--n==0){current=first=tag=null;}
	}
	public void delete(Item item) {
		if(n==0) return;
		boolean found=false;
		List list=first;
		do {
			if(list.item==item) { found=true; break; }
			list=list.next;
		}	while (list!=first);
		if(!found) { 
			System.out.println("CAN'T DELETE FROM LIST: ITEM NOT FOUND");
			return;
		}
		list.item=null;
		if(tag==list)tag=null;
		if(list==current)current=current.next;
		if(list==first)first=first.next;
		list.prev.next=list.next;
		list.next.prev=list.prev;
		list=null;
		if(--n==0){current=first=tag=null;}
	}
	public void goFirst() {
		current=first;
	} 
	public void goLast() {
		current=first.prev;
	} 
	public void goNext() {
		if(current!=null) current=current.next;
//System.out.println("goNext: current="+current);///DDD
	} 
	public void goPrev() {
		if(current!=null) current=current.prev;
	} 
	public boolean isFirst() {
		return current==first?true:false;
	}
	public boolean isLast() {
		return current==first.prev?true:false;
	}
	public void goTag() {
		current=tag;
	}
	public void tag() {
		tag=current;
	}
	public void untag() {
		tag=null;
	}
	public void goTag(int itag) {
		if(itag<1||itag>=ntags) {
			System.out.println("CAN'T GO TO TAG NUMBER "+itag);
			System.exit(1);
		}
		current=tag1;//[itag];
	}
	public void tag(int itag) {
		if(itag<1||itag>=ntags) {
			System.out.println("CAN'T USE TAG NUMBER "+itag);
			System.exit(1);
		}
		//? tag[itag]=current;
		tag1=current;
	}
	public boolean isTag(int itag) { 
		if(itag<1||itag>=ntags) {
			System.out.println("USE TAG NUMBER "+itag+" EXCEEDS MAX "+ntags);
			System.exit(1);
		}
		//? return current==tag[itag]?true:false;
		return current==tag1?true:false;
	}
	public Item getItem() {
		if(current==null)return null;
		return current.item;
	}
	public void replaceItem(Item newitem) {
		if(current==null)return;
		current.item=newitem;
	}
	public Item getTaggedItem() {
		if(tag==null)return null;
		return tag.item;
	}
	public void clean() {
		//while (number()>0) delete();
		current=first=tag=tag1=null;
		n=0;
	}
}

